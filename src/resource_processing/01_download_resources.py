from src.common.io import download_from_web, unpack_file, delete_file
from src.common.paths import resolve_dir_from_root


# GO
GO_BASIC_JSON_URL: str = "https://purl.obolibrary.org/obo/go/go-basic.json"
RGD_GAF_URL: str = "https://current.geneontology.org/annotations/rgd.gaf.gz"
MGI_GAF_URL: str = "https://current.geneontology.org/annotations/mgi.gaf.gz"

# Reactome
REACTOME_GMT_URL = "https://reactome.org/download/current/ReactomePathways.gmt.zip"
REACTOME_PATHWAYS_URL = "https://reactome.org/download/current/ReactomePathways.txt"
REACTOME_RELATIONS_URL = "https://reactome.org/download/current/ReactomePathwaysRelation.txt"


def download_go_resources():


    output_dir = resolve_dir_from_root("resources/go")

    download_from_web(GO_BASIC_JSON_URL, output_dir / "go-basic.json")

    download_from_web(RGD_GAF_URL, output_dir / "rgd.gaf.gz")

    download_from_web(MGI_GAF_URL, output_dir / "mgi.gaf.gz")


def download_go_reference_data():
    output_dir = resolve_dir_from_root("resources/reactome")

    download_from_web(REACTOME_GMT_URL, output_dir / "reactome_pathways.gmt.zip")
    unpack_file(output_dir / "reactome_pathways.gmt.zip", output_dir)
    #delete_file(output_dir / "reactome_pathways.gmt.zip")

    download_from_web(REACTOME_PATHWAYS_URL, output_dir / "reactome_pathways.txt")
    download_from_web(REACTOME_RELATIONS_URL, output_dir / "reactome_pathways_relation.txt")

