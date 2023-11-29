#! /usr/bin/env python

import argparse
import sys
import logging
import logging.config
from cellmaps_utils import logutils
from cellmaps_utils import constants
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
                        help='Minimum term size to consider for enrichment')
    parser.add_argument('--corum', default='764f7471-9b79-11ed-9a1f-005056ae23aa',
                        help='UUID for CORUM network')
    parser.add_argument('--go_cc', default='f484e8ee-0b0f-11ee-aa50-005056ae23aa',
                        help='UUID for GO-CC network')
    parser.add_argument('--hpa', default='a6a88e2d-9c0f-11ed-9a1f-005056ae23aa',
                        help='UUID for HPA network')
    parser.add_argument('--ndex_server', default='http://www.ndexbio.org',
                        help='NDEx server to use')
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
                                           corum=theargs.corum,
                                           go_cc=theargs.go_cc,
                                           hpa=theargs.hpa,
                                           ndex_server=theargs.ndex_server,
                                           hierarchy_dir=theargs.hierarchy_dir,
                                           skip_logging=theargs.skip_logging,
                                           input_data_dict=theargs.__dict__).run()
    except Exception as e:
        logger.exception('Caught exception: ' + str(e))
        return 2
    finally:
        logging.shutdown()


if __name__ == '__main__':  # pragma: no cover
    sys.exit(main(sys.argv))
