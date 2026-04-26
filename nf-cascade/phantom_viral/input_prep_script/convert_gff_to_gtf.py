import gffutils
import argparse

def convertGFFtoGTF(args):
    input_gff3 = args.input_gff
    output_gtf = args.output_gtf

    # Create a database in memory
    db = gffutils.create_db(input_gff3, dbfn=":memory:", force=True,
                            keep_order=True, merge_strategy="merge",
                            sort_attribute_values=True)

    with open(output_gtf, "w") as out:
        for gene in db.features_of_type("gene", order_by="start"):
            gene_id = gene.attributes.get("Dbxref", [gene.id])[0].replace("GeneID:", "")
            gene_name = gene.attributes.get("Name", [gene_id])[0]
            gene_type = gene.attributes.get("gene_biotype", ["protein_coding"])[0]

            # Gene line
            out.write(f"{gene.seqid}\tRefSeq\tgene\t{gene.start}\t{gene.end}\t.\t{gene.strand}\t.\t"
                      f'gene_id "{gene_id}"; gene_name "{gene_name}"; gene_type "{gene_type}";\n')

            # Look for transcripts (mRNA) under this gene
            mrnas = list(db.children(gene, featuretype="mRNA", order_by="start"))

            # Case 1: gene → mRNA → CDS
            if mrnas:
                for mrna in mrnas:
                    transcript_id = mrna.id

                    cds_features = list(db.children(mrna, featuretype="CDS", order_by="start"))
                    protein_id = cds_features[0].attributes.get("protein_id", [""])[0] if cds_features else ""

                    # Transcript line
                    out.write(f"{gene.seqid}\tRefSeq\ttranscript\t{mrna.start}\t{mrna.end}\t.\t{gene.strand}\t.\t"
                              f'gene_id "{gene_id}"; gene_name "{gene_name}"; gene_type "{gene_type}"; '
                              f'transcript_id "{transcript_id}"; protein_id "{protein_id}";\n')

                    exon_counter = 1
                    for cds in cds_features:
                        exon_id = f"exon-{transcript_id}-{exon_counter}"

                        # Exon
                        out.write(f"{gene.seqid}\tRefSeq\texon\t{cds.start}\t{cds.end}\t.\t{gene.strand}\t.\t"
                                  f'gene_id "{gene_id}"; gene_name "{gene_name}"; gene_type "{gene_type}"; '
                                  f'transcript_id "{transcript_id}"; exon_id "{exon_id}"; protein_id "{protein_id}";\n')

                        # CDS
                        frame = cds.frame if hasattr(cds, "frame") else 0
                        out.write(f"{gene.seqid}\tRefSeq\tCDS\t{cds.start}\t{cds.end}\t.\t{gene.strand}\t{frame}\t"
                                  f'gene_id "{gene_id}"; gene_name "{gene_name}"; gene_type "{gene_type}"; '
                                  f'transcript_id "{transcript_id}"; exon_id "{exon_id}"; protein_id "{protein_id}";\n')

                        exon_counter += 1

                    # Start/stop codon
                    if cds_features:
                        if gene.strand == "+":
                            start, stop = cds_features[0], cds_features[-1]
                            out.write(f"{gene.seqid}\tRefSeq\tstart_codon\t{start.start}\t{start.start+2}\t.\t+\t0\t"
                                      f'gene_id "{gene_id}"; gene_name "{gene_name}"; transcript_id "{transcript_id}"; protein_id "{protein_id}";\n')
                            out.write(f"{gene.seqid}\tRefSeq\tstop_codon\t{stop.end-2}\t{stop.end}\t.\t+\t0\t"
                                      f'gene_id "{gene_id}"; gene_name "{gene_name}"; transcript_id "{transcript_id}"; protein_id "{protein_id}";\n')
                        else:
                            start, stop = cds_features[-1], cds_features[0]
                            out.write(f"{gene.seqid}\tRefSeq\tstart_codon\t{start.end-2}\t{start.end}\t.\t-\t0\t"
                                      f'gene_id "{gene_id}"; gene_name "{gene_name}"; transcript_id "{transcript_id}"; protein_id "{protein_id}";\n')
                            out.write(f"{gene.seqid}\tRefSeq\tstop_codon\t{stop.start}\t{stop.start+2}\t.\t-\t0\t"
                                      f'gene_id "{gene_id}"; gene_name "{gene_name}"; transcript_id "{transcript_id}"; protein_id "{protein_id}";\n')

            # Case 2: gene → CDS (no mRNA)
            else:
                transcript_id = gene_id
                cds_features = list(db.children(gene, featuretype="CDS", order_by="start"))
                protein_id = cds_features[0].attributes.get("protein_id", [""])[0] if cds_features else ""

                # Transcript line (gene span)
                out.write(f"{gene.seqid}\tRefSeq\ttranscript\t{gene.start}\t{gene.end}\t.\t{gene.strand}\t.\t"
                          f'gene_id "{gene_id}"; gene_name "{gene_name}"; gene_type "{gene_type}"; '
                          f'transcript_id "{transcript_id}"; protein_id "{protein_id}";\n')

                exon_counter = 1
                if cds_features:
                    for cds in cds_features:
                        exon_id = f"exon-{transcript_id}-{exon_counter}"

                        # Exon
                        out.write(f"{gene.seqid}\tRefSeq\texon\t{cds.start}\t{cds.end}\t.\t{gene.strand}\t.\t"
                                  f'gene_id "{gene_id}"; gene_name "{gene_name}"; gene_type "{gene_type}"; '
                                  f'transcript_id "{transcript_id}"; exon_id "{exon_id}"; protein_id "{protein_id}";\n')

                        # CDS
                        frame = cds.frame if hasattr(cds, "frame") else 0
                        out.write(f"{gene.seqid}\tRefSeq\tCDS\t{cds.start}\t{cds.end}\t.\t{gene.strand}\t{frame}\t"
                                  f'gene_id "{gene_id}"; gene_name "{gene_name}"; gene_type "{gene_type}"; '
                                  f'transcript_id "{transcript_id}"; exon_id "{exon_id}"; protein_id "{protein_id}";\n')

                        exon_counter += 1

                    # Start/stop codon
                    if gene.strand == "+":
                        start, stop = cds_features[0], cds_features[-1]
                        out.write(f"{gene.seqid}\tRefSeq\tstart_codon\t{start.start}\t{start.start+2}\t.\t+\t0\t"
                                  f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; protein_id "{protein_id}";\n')
                        out.write(f"{gene.seqid}\tRefSeq\tstop_codon\t{stop.end-2}\t{stop.end}\t.\t+\t0\t"
                                  f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; protein_id "{protein_id}";\n')
                    else:
                        start, stop = cds_features[-1], cds_features[0]
                        out.write(f"{gene.seqid}\tRefSeq\tstart_codon\t{start.end-2}\t{start.end}\t.\t-\t0\t"
                                  f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; protein_id "{protein_id}";\n')
                        out.write(f"{gene.seqid}\tRefSeq\tstop_codon\t{stop.start}\t{stop.start+2}\t.\t-\t0\t"
                                  f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; protein_id "{protein_id}";\n')
                else:
                    # fallback: one exon for whole gene
                    exon_id = f"exon-{transcript_id}-1"
                    out.write(f"{gene.seqid}\tRefSeq\texon\t{gene.start}\t{gene.end}\t.\t{gene.strand}\t.\t"
                              f'gene_id "{gene_id}"; transcript_id "{transcript_id}"; exon_id "{exon_id}";\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert GFF3 to GTF (virus/human compatible)")
    parser.add_argument("input_gff", type=str, help="<Required> input gff path")
    parser.add_argument("output_gtf", type=str, help="<Required> output gtf path")
    args = parser.parse_args()
    convertGFFtoGTF(args)
