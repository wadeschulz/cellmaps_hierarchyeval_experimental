#! /usr/bin/env python

import argparse
import sys
import logging
import logging.config
from cellmaps_utils import logutils
from cellmaps_utils import constants
from cellmaps_utils.provenance import ProvenanceUtil
import cellmaps_hierarchyeval
from cellmaps_hierarchyeval.runner import CellmapshierarchyevalRunner

logger = logging.getLogger(__name__)


HIERARCHYDIR = '--hierarchy_dir'


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
    parser.add_argument('--max_fdr', type=float, default='0.05',
                        help='Maximum false discovery rate')
    parser.add_argument('--min_jaccard_index', type=float, default=0.1,
                        help='Minimum jaccard index')
    parser.add_argument('--min_comp_size', type=int, default=4,
                        help='Minimum enrichment overlap comparison size')
    parser.add_argument('--logconf', default=None,
                        help='Path to python logging configuration file in '
                             'this format: https://docs.python.org/3/library/'
                             'logging.config.html#logging-config-fileformat '
                             'Setting this overrides -v parameter which uses '
                             ' default logger. (default None)')
    parser.add_argument('--verbose', '-v', action='count', default=0,
                        help='Increases verbosity of logger to standard '
                             'error for log messages in this module. Messages are '
                             'output at these python logging levels '
                             '-v = ERROR, -vv = WARNING, -vvv = INFO, '
                             '-vvvv = DEBUG, -vvvvv = NOTSET (default no '
                             'logging)')
    parser.add_argument('--version', action='version',
                        version=('%(prog)s ' +
                                 cellmaps_hierarchyeval.__version__))

    return parser.parse_args(args)


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
    Takes a HiDeF {hierarchy_file} file from {hierarchy_dir} and runs enrichment tests for GO, CORUM, and HPA terms. 

    """.format(version=cellmaps_hierarchyeval.__version__,
               hierarchy_file=constants.HIERARCHY_NETWORK_PREFIX,
               hierarchy_dir=HIERARCHYDIR)
    
    theargs = _parse_arguments(desc, args[1:])
    theargs.program = args[0]
    theargs.version = cellmaps_hierarchyeval.__version__
    try:
        logutils.setup_cmd_logging(theargs)  
        return CellmapshierarchyevalRunner(outdir=theargs.outdir,
                                           max_fdr=theargs.max_fdr,
                                           min_jaccard_index=theargs.min_jaccard_index,
                                           min_comp_size=theargs.min_comp_size,
                                           hierarchy_dir=theargs.hierarchy_dir,
                                           input_data_dict=theargs.__dict__).run()
    except Exception as e:
        logger.exception('Caught exception: ' + str(e))
        return 2
    finally:
        logging.shutdown()


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))
