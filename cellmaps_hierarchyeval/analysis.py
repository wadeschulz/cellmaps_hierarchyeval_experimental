
import os
import re
import subprocess
import logging
from cellmaps_hierarchyeval.exceptions import CellmapshierarchyevalError
from ndex2.cx2 import CX2Network

logger = logging.getLogger(__name__)


class Hierarchy(CX2Network):
    """
    Represents an assembly of proteins in a Hierarchy
    """
    def __init__(self):
        """
        Constructor
        """
        pass


class Assembly(object):
    """
    Represents genes
    """
    def __init__(self):
        """
        Constructor
        """
        pass


class GenesetAgent(object):
    """
    Represents a Gene set analysis agent
    whose job is to consume a list of gene names
    and return a term name, confidence score, and
    analysis
    """

    NAME = 'NAME'
    """
    Term name analysis return value
    """

    CONFIDENCE = 'CONFIDENCE'
    """
    Confidence analysis return value
    """

    ANALYSIS = 'ANALYSIS'
    """
    Analysis return value
    """

    ALL = [NAME, CONFIDENCE, ANALYSIS]
    """
    List of all analysis return values
    """

    GENE_SET_TOKEN = '@@GENE_SET@@'

    def __init__(self):
        """
        Constructor
        """
        pass

    def annotate_gene_set(self, gene_names=None,
                          return_values=None):
        """
        Should be implemented by subclasses

        :param gene_names: gene symbols
        :type gene_names:
        :param return_values: Desired analysis to be returned, can be a list
                              of one or more of
        :type return_values: list
        :return:
        :rtype: dict
        """
        raise NotImplementedError('Subclasses should implement')


class OllamaCommandLineGeneSetAgent(GenesetAgent):
    """
    Runs
    """

    DEFAULT_PROMPT_FILE = 'default_prompt.txt'

    def __init__(self, prompt=None, model='llama2:latest',
                 ollama_binary='ollama'):
        """
        Constructor

        :param prompt: Prompt to pass to LLM put @@GENE_SET@@
                       into prompt to denote where gene set
                       should be inserted. If ``None`` default
                       internal prompt is used
        :type prompt: str
        """

        super().__init__()
        if prompt is None:
            logger.debug('Using default prompt')
            self._prompt = self.get_default_prompt()
        else:
            self._prompt = prompt
        self._model = model
        self._ollama_binary = ollama_binary

    def get_prompt(self):
        """
        Gets prompt used by this agent
        :return:
        """
        return self._prompt

    def get_default_prompt(self):
        """
        Gets default prompt stored with this package

        :return:
        :rtype: str
        """
        with open(os.path.join(os.path.dirname(__file__),
                  OllamaCommandLineGeneSetAgent.DEFAULT_PROMPT_FILE),
                  'r') as f:
            return f.read()

    def _run_cmd(self, cmd, cwd=None, timeout=360):
        """
        Runs command as a command line process

        :param cmd: command to run
        :type cmd: list
        :param cwd: current working directory
        :type cwd: str
        :param timeout: timeout in seconds before killing process
        :type timeout: int or float
        :raises CellMapsProvenanceError: If **raise_on_error** passed
                                         into constructor is ``True`` and
                                         process times out before completing
        :return: (return code, standard out, standard error)
        :rtype: tuple
        """
        logger.debug('Running command under ' + str(cwd) +
                     ' path: ' + str(cmd))
        p = subprocess.Popen(cmd, cwd=cwd,
                             text=True,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        try:
            out, err = p.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            logger.warning('Timeout reached. Killing process')
            p.kill()
            out, err = p.communicate()
            raise CellmapshierarchyevalError('Process timed out. '
                                             'exit code: ' +
                                             str(p.returncode) +
                                             ' stdout: ' + str(out) +
                                             ' stderr: ' + str(err))

        # Removing ending new line if value is not None
        if out is not None:
            out = out.rstrip()
        return p.returncode, out, err

    def _update_prompt_with_gene_set(self, gene_names=None):
        """
        Updates prompt inserting gene names
        :param gene_names:
        :type gene_names: list
        :return: prompt with gene names inserted
        :rtype: str
        """

        return re.sub(GenesetAgent.GENE_SET_TOKEN,
                      ','.join(gene_names), self._prompt)

    def annotate_gene_set(self, gene_names=None,
                          return_values=None):
        """
        Using prompt passed in via constructor, this call
        invokes the LLM specified by **model** set in constructor

        :param gene_names:
        :type gene_names: list
        :param return_values:
        :type return_values: list
        :return: responses from agent
        :rtype: dict
        """
        if return_values is None:
            ret_vals = GenesetAgent.ALL
        else:
            ret_vals = return_values

        updated_prompt = self._update_prompt_with_gene_set(gene_names=gene_names)

        e_code, out, err = self._run_cmd([self._ollama_binary, 'run', self._model, updated_prompt])

        return e_code, out, err

