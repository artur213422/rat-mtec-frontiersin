from typing import Any

from src.reference_processing.internal.enums import to_namespace
from src.reference_processing.internal.models import GoPathway

GO_PREFIX: str = "http://purl.obolibrary.org/obo/"
GO_NAMESPACE_PREDICATE: str = "http://www.geneontology.org/formats/oboInOwl#hasOBONamespace"



def uri_to_go_id(uri: str) -> str:
    if not uri.startswith(GO_PREFIX):
        raise ValueError(f"Expected GO URI, got: {uri}")

    value = uri.removeprefix(GO_PREFIX)

    # GO_0019882 -> GO:0019882
    return value.replace("_", ":")


def get_namespace(node: dict[str, Any]) -> str | None:
    meta = node.get("meta")

    if not isinstance(meta, dict):
        return None

    basic_property_values = meta.get("basicPropertyValues")

    if not isinstance(basic_property_values, list):
        return None

    for item in basic_property_values:
        if not isinstance(item, dict):
            continue

        pred = item.get("pred")
        val = item.get("val")

        if pred == GO_NAMESPACE_PREDICATE and isinstance(val, str):
            return val

    return None

def map_go_pathways(data: dict[str, Any]) -> list[GoPathway]:
    graphs = data.get("graphs")

    if not isinstance(graphs, list) or not graphs:
        raise ValueError("Expected JSON to contain non-empty 'graphs' list.")

    graph = graphs[0]

    if not isinstance(graph, dict):
        raise TypeError("Expected first graph to be a dictionary.")

    nodes = graph.get("nodes")

    if not isinstance(nodes, list):
        raise ValueError("Expected graph to contain 'nodes' list.")

    pathways: list[GoPathway] = []

    for node in nodes:
        if not isinstance(node, dict):
            continue

        uri = node.get("id")
        name = node.get("lbl")

        if not isinstance(uri, str):
            continue

        if not uri.startswith(GO_PREFIX + "GO_"):
            continue

        if not isinstance(name, str):
            continue

        meta = node.get("meta", {})

        # Skip obsolete terms
        if isinstance(meta, dict) and meta.get("deprecated") is True:
            continue

        namespace = get_namespace(node)

        if namespace is None:
            continue

        pathways.append(
            GoPathway(
                namespace=to_namespace(namespace),
                go_term=uri_to_go_id(uri),
                name=name,
            )
        )

    return pathways