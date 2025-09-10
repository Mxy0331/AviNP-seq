import subprocess
import argparse
import os
import pandas as pd
import sys

def run_command(command):
    """Run a shell command and check for errors."""
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        print(f"Error running command: {command}\n{stderr.decode()}")
        exit(1)
    else:
        print(f"Successfully ran command: {command}\n{stdout.decode()}")

def blat(reference_genome, input_fasta, output_psl):
    """Align the FASTA file to the reference genome using BLAT with no header in the output PSL file."""
    command = f"blat {reference_genome} {input_fasta} {output_psl} -noHead"
    run_command(command)

def extract_unaligned_regions_and_generate_report(psl_file, input_fasta, output_unaligned_fasta, report_file):
    """Extract unaligned regions from BLAT PSL output and generate a detailed report."""
    aligned_positions = {}
    best_hits = {}
    with open(psl_file, 'r') as psl:
        for line in psl:
            fields = line.strip().split('\t')
            if len(fields) < 21:
                print(f"Skipping line due to insufficient fields: {line.strip()}")
                continue
            score = int(fields[0])
            start = fields[11]
            end = fields[12]
            qsize = fields[10]
            identity = fields[0]
            chrom = fields[13]
            strand = fields[8]
            target_start = fields[15]
            target_end = fields[16]
            span = int(end) - int(start)
            query = fields[9]
            if query not in best_hits or score > best_hits[query][0]:
                best_hits[query] = [score, start, end, qsize, identity, chrom, strand, target_start, target_end, span]
            q_start = int(fields[11])
            q_end = int(fields[12])
            if query not in aligned_positions:
                aligned_positions[query] = []
            aligned_positions[query].append((q_start, q_end))
    
    report_data = [[query] + best_hits[query] for query in best_hits]
    
    unaligned_sequences = {}
    all_queries = set()
    with open(input_fasta, 'r') as fasta:
        sequence = ""
        query = ""
        for line in fasta:
            if line.startswith('>'):
                if query:
                    all_queries.add(query)
                    unaligned_sequences[query] = get_unaligned_parts(sequence, aligned_positions.get(query, []))
                query = line.strip()[1:]
                sequence = ""
            else:
                sequence += line.strip()
        if query:
            all_queries.add(query)
            unaligned_sequences[query] = get_unaligned_parts(sequence, aligned_positions.get(query, []))
    
    with open(output_unaligned_fasta, 'w') as output_fasta:
        for query, parts in unaligned_sequences.items():
            for i, part in enumerate(parts):
                output_fasta.write(f">{query}_unaligned_{i}\n{part}\n")

    # Identify sequences with no alignment
    unaligned_only = all_queries - set(best_hits.keys())
    
    # Save report as Excel file
    df = pd.DataFrame(report_data, columns=["ID", "SCORE", "START", "END", "QSIZE", "IDENTITY", "CHROM", "STRAND", "TARGET_START", "TARGET_END", "SPAN"])
    df.to_excel(report_file, index=False)
    
    # Save unaligned-only report as Excel file
    unaligned_df = pd.DataFrame(list(unaligned_only), columns=["ID"])
    unaligned_df.to_excel("unaligned_only_report.xlsx", index=False)

def get_unaligned_parts(sequence, aligned_positions):
    """Get unaligned parts of a sequence given the aligned positions."""
    if aligned_positions:
        aligned_positions.sort()
    else:
        return [sequence]  # If no aligned positions, the entire sequence is unaligned
    unaligned_parts = []
    prev_end = 0
    for start, end in aligned_positions:
        if start > prev_end:
            unaligned_parts.append(sequence[prev_end:start])
        prev_end = end
    if prev_end < len(sequence):
        unaligned_parts.append(sequence[prev_end:])
    return unaligned_parts

def main():
    reference_genome = sys.argv[1]
    input_fastq = "nomatch.fastq"
    input_fasta = "nomatch.fasta"
    output_psl = "aligned_reads.psl"
    output_unaligned_fasta = "unaligned_reads.fasta"
    report_file = "detailed_report.xlsx"

    command = f"seqkit fq2fa {input_fastq} -o {input_fasta}"
    run_command(command)

    # Step 1: Align the FASTA file to the reference genome using BLAT with no header
    print("Aligning FASTA file to reference genome using BLAT...")
    blat(reference_genome, input_fasta, output_psl)

    # Step 2: Extract unaligned regions and generate detailed report
    print("Extracting unaligned regions and generating detailed report...")
    extract_unaligned_regions_and_generate_report(output_psl, input_fasta, output_unaligned_fasta, report_file)

    # print("Pipeline completed successfully. Unaligned regions have been extracted and saved.")
    # print(f"Detailed report generated: {report_file}")
    # print("Unaligned only report generated: unaligned_only_report.xlsx")

if __name__ == "__main__":
    main()
