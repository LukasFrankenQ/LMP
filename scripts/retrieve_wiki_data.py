# -*- coding: utf-8 -*-
# Copyright 2024-2024 Lukas Franken (University of Edinburgh, Octopus Energy)
# SPDX-FileCopyrightText: : 2024-2024 Lukas Franken
#
# SPDX-License-Identifier: MIT
"""
This rule downloads data on Balancing Mechanism Units from https://query.wikidata.org

**Outputs**

- ``RESOURCES/wiki_data.csv``: collected data

"""

import logging
import pandas as pd

from SPARQLWrapper import SPARQLWrapper, JSON

from _helpers import configure_logging

logger = logging.getLogger(__name__)

if __name__ == "__main__":

    configure_logging(snakemake)

    logger.info("Retrieving BMU from wikidata.")

    sparql_url = "https://query.wikidata.org/sparql"

    # Create a SPARQLWrapper object with the defined URL
    sparql = SPARQLWrapper(sparql_url, agent='example-UA (https://example.com/; mail@example.com)')

    # Define your SPARQL query
    sparql_query = """
    SELECT DISTINCT ?item ?itemLabel ?bmrs_id ?repd_id ?lat ?lon ?capacity ?typeLabel ?instance WHERE {
    ?item wdt:P11610 ?bmrs_id.
    OPTIONAL { ?item wdt:P9891 ?repd_id. }
    OPTIONAL { ?item wdt:P2109 ?capacity. }
    OPTIONAL {
        ?item p:P625 ?point .
        ?point psv:P625 ?point_value .
        ?point_value wikibase:geoLatitude ?lat.
        ?point_value wikibase:geoLongitude ?lon.
    }
    OPTIONAL {
        ?item wdt:P31 ?type.
        ?type wdt:P279+ wd:Q159719.
    }
    OPTIONAL { ?item wdt:P31 ?instance. }
    SERVICE wikibase:label { bd:serviceParam wikibase:language "[AUTO_LANGUAGE]". }
    }
    """

    sparql.setQuery(sparql_query)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    get_quantity = lambda key, _type: [
        _type(result.get(key, {"value": "0"})["value"])
        for result in results["results"]["bindings"]
        ]
    
    wiki_df = pd.DataFrame()

    for key, _type in [
        ("lat", float),
        ("lon", float),
        ("bmrs_id", str),
        ("repd_id", str),
        ("capacity", float),
        ("item", str),
        ("itemLabel", str),
        ("typeLabel", str),
        ("instance", str),
    ]:
        wiki_df[key] = get_quantity(key, _type)

    logger.info(f"Retrieved {len(wiki_df)} BMUs from wikidata.")
    logger.info(f"Writing to {snakemake.output['wiki_data']}.")

    wiki_df.to_csv(snakemake.output["wiki_data"], index=False)