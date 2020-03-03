import utils
import TSO500
from annotate import get_annotated_snv
from reportable_variants import is_tso_snv_reportable
from vep37 import get_key_from_variant


def main():
    client = utils.get_mongo_client()
    db = utils.get_database(client,'omni')
    path = 'data/TSO500_UniqueVariants_runs1-6.csv'
    # path = 'data/TSO500_UniqueVariants_4.csv'
    var_list = TSO500.read_tso_unique_variants(path)
    outF = open("data/decisions.tsv", "w")
    h = "gene\tcdot\tpdot\tgene_category\tmutation_type\treport_status\treasons\tlack_of_reasons\t" \
        "is_protein_altering\tin_clinvar\tclinvar_significance\t" \
        "is_gain_of_function\tis_loss_of_function\thotspots\tpredicted_deleterious\tis_truncating_variants\tis_near_GOF_LOF_mutation\t" \
        "omni_gene\tomni_cdot\tomni_pdot"
    # print(h)
    outF.write(h)
    outF.write("\n")

    for index, variant in enumerate(var_list):
        key = get_key_from_variant(variant)
        annotated = get_annotated_snv(key,db)
        if not 'reasons' in annotated:
            annotated['reasons'] = []
        if not 'lack_of_reasons' in annotated:
            annotated['lack_of_reasons'] = []

        is_tso_snv_reportable(annotated)
        s = f"{annotated['gene']}\t{annotated['cdot']}\t{annotated['pdot']}\t{annotated['gene_category']}\t{annotated['mutation_type']}\t" \
            f"{annotated['report_status']}\t{annotated['reasons']}\t{annotated['lack_of_reasons']}\t" \
            f"{annotated['is_protein_altering']}\t{annotated['in_clinvar']}\t{annotated['clinvar_significance']}\t" \
            f"{annotated['is_gain_of_function']}\t{annotated['is_loss_of_function']}\t{annotated['hotspots']}\t{annotated['predicted_deleterious']}\t" \
            f"{annotated['is_truncating_variants']}\t{annotated['is_near_GOF_LOF_mutation']}\t" \
            f"{annotated['HGNC_Symbol']}\t{annotated['omni_c_dot']}\t{annotated['omni_p_dot']}"
        # print(s)
        outF.write(s)
        outF.write("\n")
    outF.close()


if __name__ == "__main__":
    main()