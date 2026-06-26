
from gseapy import prerank, Prerank

from src.common.feratures.gsea.models import GseaInput, GseaResult


def run_gsea(input: GseaInput) -> list[GseaResult]:
    print_gsea_gene_set_diagnostics(input)

    raw: Prerank = prerank(
        rnk=input.rnk.to_pandas(),
        gene_sets=input.gene_sets_dict,
        outdir=None,
        min_size=input.min_size,
        max_size=input.max_size,
        permutation_num=input.permutation_num,
    )


    results: list[GseaResult] = []

    for row in raw.res2d.to_dict("records"):
        lead_genes = (
            row["Lead_genes"].split(";")
            if row["Lead_genes"]
            else []
        )

        results.append(
            GseaResult(
                term=row["Term"],
                es=float(row["ES"]),
                nes=float(row["NES"]),
                nominal_p_value=float(row["NOM p-val"]),
                fdr_q_value=float(row["FDR q-val"]),
                leading_edge_genes=lead_genes,
            )
        )

    return results

def print_gsea_gene_set_diagnostics(input_data: GseaInput) -> None:
    gene_col: str = input_data.rnk.columns[0]

    rnk_genes: set[str] = set(
        input_data.rnk
        .select(gene_col)
        .drop_nulls()
        .to_series()
        .to_list()
    )

    dropped_too_small: list[tuple[str, int, int]] = []
    dropped_too_large: list[tuple[str, int, int]] = []
    kept: list[tuple[str, int, int]] = []

    for term, pathway_genes in input_data.gene_sets_dict.items():
        unique_pathway_genes = set(pathway_genes)
        overlap_genes = unique_pathway_genes & rnk_genes

        original_size = len(unique_pathway_genes)
        overlap_size = len(overlap_genes)

        if overlap_size < input_data.min_size:
            dropped_too_small.append((term, original_size, overlap_size))
        elif overlap_size > input_data.max_size:
            dropped_too_large.append((term, original_size, overlap_size))
        else:
            kept.append((term, original_size, overlap_size))

    print("GSEA gene set diagnostics")
    print(f"Ranked genes: {len(rnk_genes)}")
    print(f"Input gene sets: {len(input_data.gene_sets_dict)}")
    print(f"Kept gene sets: {len(kept)}")
    print(f"Dropped because overlap < min_size={input_data.min_size}: {len(dropped_too_small)}")
    print(f"Dropped because overlap > max_size={input_data.max_size}: {len(dropped_too_large)}")


    return
    if dropped_too_small:
        print("\nDropped too small:")
        for term, original_size, overlap_size in dropped_too_small:
            print(
                f"  {term} | pathway genes={original_size}, "
                f"genes in rnk={overlap_size}"
            )

    if dropped_too_large:
        print("\nDropped too large:")
        for term, original_size, overlap_size in dropped_too_large:
            print(
                f"  {term} | pathway genes={original_size}, "
                f"genes in rnk={overlap_size}"
            )