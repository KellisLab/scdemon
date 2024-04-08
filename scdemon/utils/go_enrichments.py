#!/usr/bin/python
"""Functions for getting / managing GO enrichments for modules."""
import logging

# For GO annotations (needs connection):
from gprofiler import GProfiler


def get_goterms(obj, graph_id, attr="leiden", organism="hsapiens",
                sources=["GO:CC", "GO:BP", "GO:MF",
                         "REAC", "WP", "KEGG", "CORUM"]):
    # NOTE: Needs to be able to connect to internet to run
    # TODO: Add exception/failure mode if not able to connect
    obj.gp = GProfiler(return_dataframe=True)
    mlist = obj.get_modules(graph_id)
    obj.gpres = {}
    for ll in mlist.keys():
        testlist = mlist[ll].tolist()
        logging.info(str(ll) + ": " + " ".join(testlist))
        obj.gpres[ll] = obj.gp.profile(
            organism=organism, query=testlist, sources=sources)
        logging.info(obj.gpres[ll].head())

    return obj.gpres

# TODO: saving / returning enrichments
# TODO: plotting enrichments
# TODO: save enrichments for each graph/config
