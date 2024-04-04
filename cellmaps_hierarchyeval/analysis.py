
import os
import re
import subprocess
import random
import logging
from cellmaps_hierarchyeval.exceptions import CellmapshierarchyevalError
from ndex2.cx2 import CX2Network

logger = logging.getLogger(__name__)


class Hierarchy(object):
    """
    Represents an assembly of proteins in a Hierarchy
    """
    def __init__(self, hierarchy=None, interactome=None,
                 ndex_username=None, ndex_password=None):
        """
        Constructor
        :param hierarchy: Hierarchy
        :type hierarchy: :py:class:`~ndex2.cx2.CX2Network`
        :param interactome: Parent interactome
        :type interactome: :py:class:`~ndex2.cx2.CX2Network`
        :param ndex_username: NDEx username to use when connecting
                              to NDEx to obtain interactomes from hierarchy
        :type ndex_username: str
        :param ndex_password: NDEx password to use when connecting
                              to NDEx to obtain interactomes from hierarchy
        :type ndex_password: str
        """
        self._hierarchy = hierarchy
        self._interactome = interactome
        self._ndex_username = ndex_username
        self._ndex_password = ndex_password

    def get_next_assembly(self):
        """
        Generator that gets next assembly in hierarchy

        :return:
        :rtype: :py:class:`~cellmaps_hierarchyeval.assembly.Assembly`
        """
        raise NotImplementedError('not done yet')

        # will iterate across hierarchy CX2Network and create an assembly
        # object which contains list of gene names and needed info to link
        # back to this hierarchy node

        # for node_id, node_data in self._hierarchy.get_nodes().items():
            #yield X


class Assembly(object):
    """
    Represents assembly in a hierarchy
    """
    def __init__(self, node_id=None, gene_names=None):
        """
        Constructor

        :param node_id: Id of hierarchy node
        :type node_id: int
        :param gene_names: list of gene names
        :type gene_names: list
        """
        self._node_id = node_id
        self._gene_names = gene_names

    def get_assembly_name(self):
        """
        Gets name of assembly
        :return:
        """
        return None

    def set_assembly_name(self):
        """
        Sets assembly name
        :return:
        """
        pass

    def get_node_id(self):
        """
        Gets node id

        :return:
        """
        return self._node_id

    def get_gene_names(self):
        """
        Gets gene names

        :return:
        """
        return self._gene_names


class GenesetAgent(object):
    """
    Represents a Gene set analysis agent
    whose job is to consume a list of gene names
    and return a term name, confidence score, and
    analysis
    """
    GENE_SET_TOKEN = 'GENE_SET'

    def __init__(self, attribute_name_prefix=None):
        """
        Constructor
        """
        self._attribute_name_prefix = attribute_name_prefix

    def annotate_gene_set(self, gene_names=None):
        """
        Should be implemented by subclasses

        :param gene_names: gene symbols
        :type gene_names:
        :return:
        :rtype: tuple
        """
        raise NotImplementedError('Subclasses should implement')

    def get_attribute_name_prefix(self):
        """
        Gets suggested attribute name prefix

        :return:
        :rtype: str
        """
        return self._attribute_name_prefix


class FakeGeneSetAgent(GenesetAgent):
    """
    Fake geneset agent that generates random numbers for values
    """
    def __init__(self, random_seed=None, attribute_name_prefix=None):
        """
        Constructor
        :param random_seed:
        """
        super().__init__(attribute_name_prefix=attribute_name_prefix)
        random.seed(random_seed)
        if self._attribute_name_prefix is None:
            self._attribute_name_prefix = 'fake_' + str(random.random()) + '::'

    def annotate_gene_set(self, gene_names=None):
        """

        :param gene_names:
        :return:
        """
        return 'Fake ' +\
               str(random.randint(0, 1000)), random.random(), 'Fake full text' +\
               str(random.randint(0, 1000))


class OllamaCommandLineGeneSetAgent(GenesetAgent):
    """
    Runs
    """

    DEFAULT_PROMPT_FILE = 'default_prompt.txt'

    def __init__(self, prompt=None, model='llama2:latest',
                 ollama_binary='/usr/local/bin/ollama',
                 attribute_name_prefix=None):
        """
        Constructor

        :param prompt: Prompt to pass to LLM put @@GENE_SET@@
                       into prompt to denote where gene set
                       should be inserted. If ``None`` default
                       internal prompt is used
        :type prompt: str
        """

        super().__init__(attribute_name_prefix=attribute_name_prefix)
        if prompt is None:
            logger.debug('Using default prompt')
            self._prompt = self.get_default_prompt()
        else:
            self._prompt = prompt
        self._model = model
        self._ollama_binary = ollama_binary
        if self._attribute_name_prefix is None:
            self._attribute_name_prefix = 'ollama_' + str(self._model) + '::'

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

        return self._prompt.format(GENE_SET=','.join(gene_names))

    def annotate_gene_set(self, gene_names=None):
        """
        Using prompt passed in via constructor, this call
        invokes the LLM specified by **model** set in constructor

        :param gene_names: Genes to analyze
        :type gene_names: list
        :raises CellmapshierarchyevalError: If LLM failed to run
        :return: ('process name (score)', full output from LLM)
        :rtype: tuple
        """
        updated_prompt = self._update_prompt_with_gene_set(gene_names=gene_names)

        e_code, out, err = self._run_cmd([self._ollama_binary, 'run',
                                          self._model,
                                          updated_prompt])
        if e_code != 0:
            raise CellmapshierarchyevalError('Received non zero exit code + ' +
                                             str(e_code) +
                                             ' calling ' +
                                             str(self._ollama_binary) +
                                             '\nstdout: ' + str(out) +
                                             'stderr\n' + str(err))

        process_name = None
        confidence = None
        if out is not None:
            for line in out.split('\n'):
                if line.startswith('Process: '):
                    process_name = line[line.index(':')+2:]
                if line.startswith('Confidence Score: '):
                    confidence = line[line.index(':')+2:]
        else:
            logger.info('LLM output is None')

        return process_name, confidence, out

