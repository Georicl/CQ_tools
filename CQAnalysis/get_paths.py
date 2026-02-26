from os import path
import os

def get_paths(f_bam, m_bam, fasta_path_get, output_path_get_path) -> dict:
    """
    Centralized path management for the analysis pipeline.
    
    Args:
        f_bam (str): Path to female BAM file.
        m_bam (str): Path to male BAM file.
        fasta_path_get (str): Path to reference FASTA file.
        output_path_get_path (str): Base output directory.
        
    Returns:
        dict: A dictionary containing absolute paths for all intermediate and final files.
            - fa_path: Reference FASTA.
            - fai_path: FASTA Index (.fai).
            - chromosome_length: TSV file with chrom names and lengths.
            - windows_tsv: BED/TSV file containing genomic windows.
            - f_cov: Coverage output for female sample.
            - m_cov: Coverage output for male sample.
            - f_m_merge: Merged coverage data from both samples.
            - filtered_cq: Final filtered results based on CQ value.
            - temp_dir: Directory for intermediate temporary files.
    """
    temp_dir = path.join(output_path_get_path, "temp_dir")
    if not path.exists(temp_dir):
        os.makedirs(temp_dir)
    return {
        "fa_path": f"{fasta_path_get}",
        "log_path": path.join(output_path_get_path, "CQ.log"),
        "fai_path": f"{fasta_path_get}.fai",
        "chromosome_length": path.join(output_path_get_path, "chromosome_length.tsv"),
        "windows_tsv": path.join(output_path_get_path, "chromosome_windows.tsv"),
        "f_cov": path.join(temp_dir, "F_windows_coverage.tsv"),
        "m_cov": path.join(temp_dir, "M_windows_coverage.tsv"),
        "f_m_merge": path.join(temp_dir, "F_M_merge.tsv"),
        "filtered_cq": path.join(output_path_get_path, "F_M_CQ.filter.tsv"),
        "temp_dir": path.join(output_path_get_path, "temp_dir"),
        "f_bam": f_bam,
        "m_bam": m_bam,
    }