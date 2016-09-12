#!/usr/bin/env python2

from __future__ import print_function

import os
import sys
import tempfile

try:
    import gringo
except ImportError:
    print("!!! gringo python module is not installed. Some features will fail.",
            file=sys.stderr)

from caspots.input import *
from caspots.formats import *
from caspots.utils import *
from caspots.dataset import *
from caspots import identify
from caspots import modelchecking
from caspots.model import *

def prepare_output_file(path):
    dbg("# writing to '%s'" % path)

def debug_file(args, name):
    path = os.path.join(args.debug_dir,name)
    prepare_output_file(path)
    return path

def read_pkn(args):
    termset, graph = termset_of_sif(args.pkn)
    if args.debug:
        dfile = debug_file(args,
            "%s.lp"%os.path.basename(args.pkn).replace(".sif", ""))
        termset.to_file(dfile, pprint=True)
    return termset, graph

def dataset_name(args):
    return os.path.basename(args.dataset).replace(".csv", "")

def read_dataset(args, graph):
    termset = termset_of_dataset(args.dataset, graph, factor=args.factor)
    if args.debug:
        dfile = debug_file(args, "%s.lp"%dataset_name(args))
        termset.to_file(dfile, pprint=True)
    return termset

def read_networks(args):
    reader = component.getUtility(core.ICsvReader)
    reader.read(args.networks)
    networks = list(ILogicalNetworkList(reader))[args.range_from:]
    if args.range_length:
        networks = networks[:args.range_length]
    return networks

def read_domain(args, pkn, dataset, outf):
    if args.networks:
        networks = read_networks(args)
        out = domain_of_networks(networks, pkn, dataset)
        with open(outf, "w") as fd:
            fd.write(out)
        return outf
    else:
        return None

def do_pkn2lp(args):
    read_pkn(args)[0].to_file(args.output, pprint=True)

def do_midas2lp(args):
    _, graph = read_pkn(args)
    read_dataset(args, graph).to_file(args.output, pprint=True)

def do_results2lp(args):
    pkn, graph = read_pkn(args)
    idataset = read_dataset(args, graph)
    dataset = Dataset(dataset_name(args))
    dataset.feed_from_asp(idataset.to_str())
    networks = read_networks(args)
    out = domain_of_networks(networks, pkn, dataset)
    print(out)

def do_mse(args):
    pkn, graph = read_pkn(args)
    idataset = read_dataset(args, graph)
    termset = pkn.union(idataset)
    dataset = Dataset(dataset_name(args))
    dataset.feed_from_asp(idataset.to_str())

    fd, domainlp = tempfile.mkstemp(".lp")
    os.close(fd)
    domain = read_domain(args, pkn, dataset, domainlp)

    identifier = identify.ASPSolver(termset, args, domain=domain)


    first = True
    exact = False
    for sample in identifier.solution_samples():
        (mse0, mse) = sample.mse()
        if first:
            print("MSE_discrete = %s" % mse0)
            if mse0 == mse:
                print("MSE_sample >= MSE_discrete")
            else:
                print("MSE_sample >= %s" % mse)

        if args.check_exact:
            model = sample.model()
            trace = sample.trace(dataset)

            fd, smvfile = tempfile.mkstemp(".smv")
            os.close(fd)
            if args.debug:
                dbg("# running NuSMV on %s" % smvfile)
            exact = modelchecking.verify(trace, model, smvfile, args.semantics)
            if args.debug:
                dbg("# NuSMV says %s" % exact)
            else:
                os.unlink(smvfile)

            if exact:
                break
        else:
            break
        first = False
    if args.check_exact:
        if exact:
            print("MSE_sample is exact")
        else:
            print("MSE_sample may be under-estimated (no True Positive found)")
    os.unlink(domainlp)


def do_identify(args):
    pkn, graph = read_pkn(args)
    idataset = read_dataset(args, graph)
    termset = pkn.union(idataset)

    dataset = Dataset(dataset_name(args))
    dataset.feed_from_asp(idataset.to_str())

    fd, domainlp = tempfile.mkstemp(".lp")
    os.close(fd)
    domain = read_domain(args, pkn, dataset, domainlp)

    identifier = identify.ASPSolver(termset, args, domain=domain)

    sif = component.getUtility(core.IFileReader)
    sif.read(args.pkn)
    graph = core.IGraph(sif)
    names = component.getUtility(core.ILogicalNames)
    names.load(graph)

    networks = core.LogicalNetworkSet()

    c = {
        "found": 0,
        "tp": 0,
    }
    def show_stats(output=sys.stderr):
        if args.true_positives:
            output.write("%d solution(s) / %d true positives\r" % (c["found"], c["tp"]))
        else:
            output.write("%d solution(s)\r" % c["found"])
        output.flush()

    def is_true_positive(bid, facts):
        model = BNModel(bid)
        model.feed_from_asp(facts)
        fd, smvfile = tempfile.mkstemp(".smv")
        os.close(fd)
        exact = modelchecking.verify(dataset, model, smvfile, args.semantics)
        os.unlink(smvfile)
        return exact

    def on_model(model):
        c["found"] += 1
        facts = map(str, model.atoms(gringo.Model.SHOWN))
        skip = False
        if args.true_positives:
            if is_true_positive(c["found"], identify.parse_answer(", ".join(facts))):
                c["tp"] += 1
            else:
                skip = True
        show_stats()
        if skip:
            return
        answer = asp.AnswerSet(facts)
        network = core.ILogicalNetwork(asp.ITermSet(answer))
        networks.add(network, update_names=False)

    identifier.solutions(on_model, limit=args.limit, force_weight=args.force_weight)
    print("%d solution(s) for the over-approximation" % c["found"])
    if args.true_positives:
        print("%d/%d true positives [rate: %0.2f%%]" \
            % (c["tp"], c["found"], (100.*c["tp"])/c["found"]))

    writer = core.ICsvWriter(networks)
    writer.write(args.output)
    os.unlink(domainlp)


def do_validate(args):
    pkn, graph = read_pkn(args)
    idataset = read_dataset(args, graph)

    dataset = Dataset(dataset_name(args))
    dataset.feed_from_asp(idataset.to_str())

    reader = component.getUtility(core.ICsvReader)
    reader.read(args.networks)
    networks = list(ILogicalNetworkList(reader))[args.range_from:]
    if args.range_length:
        networks = networks[:args.range_length]
    nb = len(networks)

    def is_true_positive(bid, terms):
        model = BNModel(bid)
        model.feed_from_asp(identify.parse_answer(terms.to_str()))
        fd, smvfile = tempfile.mkstemp(".smv")
        os.close(fd)
        exact = modelchecking.verify(dataset, model, smvfile, args.semantics)
        if args.debug:
            dbg("# %s" % smvfile)
        else:
            os.unlink(smvfile)
        return exact

    tp = 0
    c = 0
    for network in networks:
        c += 1
        sys.stderr.write("%d/%d... " % (c,nb))
        sys.stderr.flush()
        if is_true_positive(c, asp.ITermSet(network)):
            tp += 1
        sys.stderr.write("%d/%d true positives\r" % (tp,c))
    res = "%d/%d true positives [rate: %0.2f%%]" % (tp, nb, (100.*tp)/nb)
    print(res)
    if args.tee:
        with open(args.tee, "w") as f:
            f.write("%s\n" % res)


from argparse import ArgumentParser

def run():

    parser = ArgumentParser(prog=sys.argv[0])
    parser.add_argument("--debug", action="store_true", default=False)
    parser.add_argument("--debug-dir", type=str, default=tempfile.gettempdir())
    subparsers = parser.add_subparsers(help="commands help")

    identify_parser = ArgumentParser(add_help=False)
    identify_parser.add_argument("--family", choices=["all", "subset", "mincard"],
                                    default="subset",
                                    help="result family (default: subset)")
    identify_parser.add_argument("--mincard-tolerance", type=int, default=0,
                                    help="consider (subset minimal) solutions with cardinality at most tolerance + the minimum cardinality")
    identify_parser.add_argument("--weight-tolerance", type=int, default=0,
                                    help="consider (subset minimal) solutions with weight at most tolerance + the minimum weight")
    identify_parser.add_argument("--enum-traces", action="store_true",
                                    default=False,
                                    help="enumerate over traces")
    identify_parser.add_argument("--fully-controllable", action="store_true",
                                    help="only consider BNs where all nodes have a stimulus in their ancestors (default)")
    identify_parser.add_argument("--no-fully-controllable", action="store_false", dest="fully_controllable",
                                    help="do not only consider BNs where all nodes have a stimulus in their ancestors")
    identify_parser.set_defaults(fully_controllable=True)
    identify_parser.add_argument("--force-weight", type=int, default=None,
                                    help="Force the maximum weight of a solution")
    identify_parser.add_argument("--force-size", type=int, default=None,
                                    help="Force the maximum size of a solution")
    modelchecking_p = ArgumentParser(add_help=False)
    modelchecking_p.add_argument("--semantics",
        choices=modelchecking.MODES, default=modelchecking.U_GENERAL,
        help="Updating mode of the Boolean network (default: %s)" \
            % modelchecking.U_GENERAL)

    pkn_parser = ArgumentParser(add_help=False)
    pkn_parser.add_argument("pkn", help="Prior knowledge network (sif format)")
    dataset_parser = ArgumentParser(add_help=False)
    dataset_parser.add_argument("dataset", help="Dataset (midas csv format)")
    dataset_parser.add_argument("--factor", type=int, default=100,
                                    help="discretization factor (default: 100)")

    networks_parser = ArgumentParser(add_help=False)
    networks_parser.add_argument("--range-from", type=int, default=0,
        help="Validate only networks from given row (starting at 0)")
    networks_parser.add_argument("--range-length", type=int, default=0,
        help="Number of networks to validate (0 means all)")

    domain_parser = ArgumentParser(add_help=False,
        parents=[networks_parser])
    domain_parser.add_argument("--networks", help="Networks to as domain (.csv)")

    parser_pkn2lp = subparsers.add_parser("pkn2lp",
        help="Export PKN (sif format) to ASP (lp format)",
        parents=[pkn_parser])
    parser_pkn2lp.add_argument("output", help="Output file (.lp format)")
    parser_pkn2lp.set_defaults(func=do_pkn2lp)

    parser_midas2lp = subparsers.add_parser("midas2lp",
        help="Export dataset (midas csv format) to ASP (lp format)",
        parents=[pkn_parser, dataset_parser])
    parser_midas2lp.add_argument("output", help="Output file (.lp format)")
    parser_midas2lp.set_defaults(func=do_midas2lp)

    parser_results2lp = subparsers.add_parser("results2lp",
        help="Export results to ASP (lp format)",
        parents=[pkn_parser, dataset_parser, networks_parser])
    parser_results2lp.add_argument("networks", help="Networks file (.csv format)")
    parser_results2lp.set_defaults(func=do_results2lp)

    parser_mse = subparsers.add_parser("mse",
        help="Compute the best MSE",
        parents=[pkn_parser, dataset_parser, identify_parser,
                    modelchecking_p, domain_parser])
    parser_mse.add_argument("--check-exact", action="store_true", default=False,
                            help="look for a true positive with the computed MSE")
    parser_mse.set_defaults(func=do_mse)

    parser_identify = subparsers.add_parser("identify",
        help="Identify all the best Boolean networks",
        parents=[pkn_parser, dataset_parser, identify_parser,
                    modelchecking_p, domain_parser])
    parser_identify.add_argument("--true-positives", default=False, action="store_true",
        help="filter solutions to keep only true positives (exact identification)")
    parser_identify.add_argument("--limit", default=0, type=int,
        help="Limit the number of solutions")
    parser_identify.add_argument("output", help="output file (csv format)")
    parser_identify.set_defaults(func=do_identify)

    parser_validate = subparsers.add_parser("validate",
        help="Compute the true positive rate of *exactly* identified networks",
        parents=[pkn_parser, dataset_parser,
                    modelchecking_p])
    parser_validate.add_argument("--range-from", type=int, default=0,
        help="Validate only networks from given row (starting at 0)")
    parser_validate.add_argument("--range-length", type=int, default=0,
        help="Number of networks to validate (0 means all)")
    parser_validate.add_argument("--tee", type=str, default=None,
        help="Output result to given file (in addition to stdout).")
    parser_validate.add_argument("networks", help="network set (csv format)")
    parser_validate.set_defaults(func=do_validate)

    args = parser.parse_args()
    args.func(args)

