#! /usr/bin/env python
import json
import os
import argparse
import sys
import logging
import logging.config
from cellmaps_utils import logutils
from cellmaps_utils import constants
import cellmaps_hierarchyeval
from cellmaps_hierarchyeval.runner import CellmapshierarchyevalRunner
from cellmaps_hierarchyeval.analysis import OllamaCommandLineGeneSetAgent
from cellmaps_hierarchyeval.analysis import OllamaRestServiceGenesetAgent
from cellmaps_hierarchyeval.analysis import FakeGeneSetAgent

logger = logging.getLogger(__name__)


HIERARCHYDIR = '--hierarchy_dir'
PATH_TO_OLLAMA = '/usr/local/bin/ollama'


def _parse_arguments(desc, args):
    """
    Parses command line arguments

    :param desc: description to display on command line
    :type desc: str
    :param args: command line arguments usually :py:func:`sys.argv[1:]`
    :type args: list
    :return: arguments parsed by :py:mod:`argparse`
    :rtype: :py:class:`argparse.Namespace`
    """
    parser = argparse.ArgumentParser(description=desc,
                                     formatter_class=constants.ArgParseFormatter)
    parser.add_argument('outdir', help='Output directory')
    parser.add_argument(HIERARCHYDIR, required=True,
                        help='Directory where hierarchy was generated')
    parser.add_argument('--max_fdr', type=float, default=CellmapshierarchyevalRunner.MAX_FDR,
                        help='Maximum false discovery rate')
    parser.add_argument('--min_jaccard_index', type=float,
                        default=CellmapshierarchyevalRunner.MIN_JACCARD_INDEX,
                        help='Minimum jaccard index')
    parser.add_argument('--min_comp_size', type=int, default=CellmapshierarchyevalRunner.MIN_COMP_SIZE,
                        help='Minimum term size to consider for enrichment')
    parser.add_argument('--corum', default=CellmapshierarchyevalRunner.CORUM,
                        help='UUID for CORUM network')
    parser.add_argument('--go_cc', default=CellmapshierarchyevalRunner.GO_CC,
                        help='UUID for GO-CC network')
    parser.add_argument('--hpa', default=CellmapshierarchyevalRunner.HPA,
                        help='UUID for HPA network')
    parser.add_argument('--ndex_server', default=CellmapshierarchyevalRunner.NDEX_SERVER,
                        help='NDEx server to use')
    parser.add_argument('--skip_term_enrichment', action='store_true',
                        help='If set, SKIP enrichment against networks set '
                             'via --corum, --go_cc, --hpa')
    parser.add_argument('--ollama', default=PATH_TO_OLLAMA,
                        help='Path to ollama command line binary or REST service. '
                             'If value starts with http it is assumed to be a REST '
                             'url and all prompts will be passed to service. For'
                             'REST url the suffix api/generate must be appended. '
                             'Example: http://foo/api/generate '
                             'NOTE: ollama integration with this tool is '
                             'EXPERIMENTAL and interface may be '
                             'changed or removed in the future ')
    parser.add_argument('--ollama_user',
                        help='Username to pass as basic auth to ollama REST '
                             'service')
    parser.add_argument('--ollama_password',
                        help='Password to pass via basic autho to ollama REST '
                             'service')
    parser.add_argument('--ollama_prompts', nargs='+',
                        help='Comma delimited value of format <MODEL NAME> or '
                             '<MODEL NAME>,<PROMPT> '
                             'where <PROMPT> can be path to prompt file or prompt to '
                             'run. For insertion of gene set please include {GENE_SET} '
                             'in prompt and tell LLM to put Process: <name> on first line '
                             'with name assigned to assembly and Confidence Score: <score> '
                             'on 2nd line with confidence in the name given. '
                             'If just <MODEL NAME> is set, then default prompt is used with '
                             'model specified. '
                             'NOTE: if <MODEL NAME> is set to FAKE then a completely fake '
                             ' agent will be used. Also note: ollama integration with this '
                             'tool is EXPERIMENTAL and interface may be '
                             'changed or removed in the future ')
    parser.add_argument('--provenance',
                        help='Path to file containing provenance '
                             'information about input files in JSON format. '
                             'This is required if inputdir does not contain '
                             'ro-crate-metadata.json file.')
    parser.add_argument('--name',
                        help='Name of this run, needed for FAIRSCAPE. If '
                             'unset, name value from specified '
                             'by --hierarchy_dir directory or provenance file will be used')
    parser.add_argument('--organization_name',
                        help='Name of organization running this tool, needed '
                             'for FAIRSCAPE. If unset, organization name specified '
                             'in --hierarchy_dir directory or provenance file will be used')
    parser.add_argument('--project_name',
                        help='Name of project running this tool, needed for '
                             'FAIRSCAPE. If unset, project name specified '
                             'in --hierarchy_dir directory or provenance file will be used')
    parser.add_argument('--skip_logging', action='store_true',
                        help='If set, output.log, error.log '
                             'files will not be created')
    parser.add_argument('--logconf', default=None,
                        help='Path to python logging configuration file in '
                             'this format: https://docs.python.org/3/library/'
                             'logging.config.html#logging-config-fileformat '
                             'Setting this overrides -v parameter which uses '
                             ' default logger. (default None)')
    parser.add_argument('--verbose', '-v', action='count', default=1,
                        help='Increases verbosity of logger to standard '
                             'error for log messages in this module. Messages are '
                             'output at these python logging levels '
                             '-v = WARNING, -vv = INFO, '
                             '-vvv = DEBUG, -vvvv = NOTSET (default ERROR '
                             'logging)')
    parser.add_argument('--version', action='version',
                        version=('%(prog)s ' +
                                 cellmaps_hierarchyeval.__version__))

    return parser.parse_args(args)


def get_ollama_geneset_agents(ollama=PATH_TO_OLLAMA, ollama_prompts=None,
                              username=None, password=None):
    """
    Parses **ollama_prompts** from argparse and creates geneset agents

    :param ollama: Path to ollama binary or REST service
    :type ollama: str
    :param ollama_prompts:
    :type ollama_prompts: list
    :return:
    """
    if ollama_prompts is None:
        return None

    res = []
    use_rest_service = False
    if ollama.startswith('http'):
        logger.info('For all agents, using ollama REST service: ' +
                    str(ollama))
        if not ollama.endswith('api/generate'):
            logger.warning(str(ollama) +
                           ' does not end with api/generate and may not work.')
        use_rest_service = True

    for o_prompt in ollama_prompts:
        model, prompt = get_model_prompt_from_string(o_prompt)
        if model.lower() == 'fake':
            logger.debug('Creating FAKE geneset agent')
            res.append(FakeGeneSetAgent())
            continue

        logger.debug('Creating ollama geneset agent for model: ' + str(model))
        if use_rest_service is True:
            agent = OllamaRestServiceGenesetAgent(rest_url=ollama, username=username,
                                                  password=password,
                                                  model=model, prompt=prompt)
        else:
            agent = OllamaCommandLineGeneSetAgent(ollama_binary=ollama,
                                                  model=model, prompt=prompt)
        res.append(agent)
    return res


def get_model_prompt_from_string(o_prompt):
    """
    Given argument from --ollama_prompts flag extract
    model and prompt which can be in following formats:

    <MODEL>
    or
    <MODEL>,<PROMPT>

    Where <MODEL> will always just be a string, but <PROMPT>
    can be a string or a path to a file

    :param o_prompt: argument passed to --ollama_prompts
    :type o_prompt: str
    :return: model, prompt
    :rtype: tuple
    """
    split_prompt = o_prompt.split(',')
    model = split_prompt[0]
    prompt = None
    if len(split_prompt) > 1:
        raw_prompt = split_prompt[1]
        if os.path.isfile(raw_prompt):
            with open(raw_prompt, 'r') as f:
                prompt = f.read()
        else:
            prompt = raw_prompt

    return model, prompt


def main(args):
    """
    Main entry point for program

    :param args: arguments passed to command line usually :py:func:`sys.argv[1:]`
    :type args: list

    :return: return value of :py:meth:`cellmaps_hierarchyeval.runner.CellmapshierarchyevalRunner.run`
             or ``2`` if an exception is raised
    :rtype: int
    """
    desc = """
    Version {version}
    Takes a HiDeF {hierarchy_file} file from {hierarchy_dir} and runs
    enrichment tests for GO, CORUM, and HPA terms.

    Also includes EXPERIMENTAL support for invocation of LLMs via Ollama command
    line or Ollama REST service.

    To use see --ollama and --ollama_prompts flags

    """.format(version=cellmaps_hierarchyeval.__version__,
               hierarchy_file=constants.HIERARCHY_NETWORK_PREFIX,
               hierarchy_dir=HIERARCHYDIR)

    theargs = _parse_arguments(desc, args[1:])
    theargs.program = args[0]
    theargs.version = cellmaps_hierarchyeval.__version__

    if theargs.provenance is not None:
        with open(theargs.provenance, 'r') as f:
            json_prov = json.load(f)
    else:
        json_prov = None

    try:
        logutils.setup_cmd_logging(theargs)

        ollama_prompts = get_ollama_geneset_agents(ollama=theargs.ollama,
                                                   ollama_prompts=theargs.ollama_prompts,
                                                   username=theargs.ollama_user,
                                                   password=theargs.ollama_password)

        return CellmapshierarchyevalRunner(outdir=theargs.outdir,
                                           max_fdr=theargs.max_fdr,
                                           min_jaccard_index=theargs.min_jaccard_index,
                                           min_comp_size=theargs.min_comp_size,
                                           corum=theargs.corum,
                                           go_cc=theargs.go_cc,
                                           hpa=theargs.hpa,
                                           ndex_server=theargs.ndex_server,
                                           geneset_agents=ollama_prompts,
                                           name=theargs.name,
                                           organization_name=theargs.organization_name,
                                           project_name=theargs.project_name,
                                           hierarchy_dir=theargs.hierarchy_dir,
                                           skip_term_enrichment=theargs.skip_term_enrichment,
                                           skip_logging=theargs.skip_logging,
                                           input_data_dict=theargs.__dict__,
                                           provenance=json_prov).run()
    except Exception as e:
        logger.exception('Caught exception: ' + str(e))
        return 2
    finally:
        logging.shutdown()


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))
