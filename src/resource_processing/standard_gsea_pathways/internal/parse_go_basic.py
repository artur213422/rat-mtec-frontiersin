from typing import Any

from src.resource_processing.standard_gsea_pathways.internal.enums import to_namespace
from src.resource_processing.standard_gsea_pathways.internal.models import GoPathway


GO_PREFIX = "http://purl.obolibrary.org/obo/"
GO_NAMESPACE_PREDICATE = "http://www.geneontology.org/formats/oboInOwl#hasOBONamespace"


# Relationships through which positive GO annotations propagate upward.
PARENT_PREDICATES = {
    "is_a",
    "http://purl.obolibrary.org/obo/BFO_0000050",  # part_of
}


def uri_to_go_id(uri: str) -> str:
    if not uri.startswith(GO_PREFIX + "GO_"):
        raise ValueError(f"Expected GO URI, got: {uri}")

    return uri.removeprefix(GO_PREFIX).replace("_", ":")


def get_namespace(node: dict[str, Any]) -> str | None:
    meta = node.get("meta")

    if not isinstance(meta, dict):
        return None

    values = meta.get("basicPropertyValues")

    if not isinstance(values, list):
        return None

    for item in values:
        if not isinstance(item, dict):
            continue

        if item.get("pred") == GO_NAMESPACE_PREDICATE:
            value = item.get("val")

            if isinstance(value, str):
                return value

    return None


def is_go_uri(value: object) -> bool:
    return isinstance(value, str) and value.startswith(GO_PREFIX + "GO_")


def map_go_pathways(data: dict[str, Any]) -> list[GoPathway]:
    graphs = data.get("graphs")

    if not isinstance(graphs, list) or not graphs:
        raise ValueError("Expected JSON to contain a non-empty 'graphs' list.")

    graph = graphs[0]

    if not isinstance(graph, dict):
        raise TypeError("Expected the first graph to be a dictionary.")

    nodes = graph.get("nodes")
    edges = graph.get("edges")

    if not isinstance(nodes, list):
        raise ValueError("Expected graph to contain a 'nodes' list.")

    if not isinstance(edges, list):
        raise ValueError("Expected graph to contain an 'edges' list.")

    pathways_by_id: dict[str, GoPathway] = {}

    # First map all valid, non-obsolete GO nodes.
    for node in nodes:
        if not isinstance(node, dict):
            continue

        uri = node.get("id")
        name = node.get("lbl")

        if not is_go_uri(uri) or not isinstance(name, str):
            continue

        meta = node.get("meta", {})

        if isinstance(meta, dict) and meta.get("deprecated") is True:
            continue

        namespace_value = get_namespace(node)

        if namespace_value is None:
            continue

        go_id = uri_to_go_id(uri)

        pathways_by_id[go_id] = GoPathway(
            namespace=to_namespace(namespace_value),
            go_term=go_id,
            name=name,
        )

    # In OBO Graphs:
    #   sub = child
    #   obj = parent
    for edge in edges:
        if not isinstance(edge, dict):
            continue

        child_uri = edge.get("sub")
        parent_uri = edge.get("obj")
        predicate = edge.get("pred")

        if predicate not in PARENT_PREDICATES:
            continue

        if not is_go_uri(child_uri) or not is_go_uri(parent_uri):
            continue

        child_id = uri_to_go_id(child_uri)
        parent_id = uri_to_go_id(parent_uri)

        child = pathways_by_id.get(child_id)
        parent = pathways_by_id.get(parent_id)

        # Ignore obsolete, missing, or external nodes.
        if child is None or parent is None:
            continue

        # GO hierarchy propagation should stay within one ontology.
        if child.namespace != parent.namespace:
            continue

        child.parent_go_terms.append(parent_id)

    # Make the output deterministic and remove duplicate edges.
    for pathway in pathways_by_id.values():
        pathway.parent_go_terms = sorted(set(pathway.parent_go_terms))

    return sorted(
        pathways_by_id.values(),
        key=lambda pathway: pathway.go_term,
    )