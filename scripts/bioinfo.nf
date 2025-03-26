import org.apache.commons.lang.RandomStringUtils

input_pdb_dir = "../pdb_processed/*.pdb"
pdb2fasta_path = "pdb2fasta"
target_fasta_dir = "../temp/out_fasta"
target_search_dir = "../temp/out_fasta_hits"
target_search_dir_msa = "../temp/out_fasta_hits_msa_input"
target_search_dir_msa_out = "../temp/out_fasta_hits_msa_output"
filter_blast_script_path = "filter_blast.py"
homo_fasta_path = "../data/UP000005640_9606.fasta"
extract_script_path = "extract.sh"
homology_count_output_path = "../outputs/homology_counts.txt"


ch1 = Channel.fromPath(input_pdb_dir)

process fasta_extractor
{
    errorStrategy 'retry'
    maxErrors 1000
    maxForks 2
    
    input:
    val(id) from ch1
    
    output:
    val(id) into ch2
    
    script:
    def base_filename = id.toString().split("/")[-1]
    def base = "mkdir -p ${target_fasta_dir} && ${pdb2fasta_path} ${id} > ${target_fasta_dir}/${base_filename}.fa"
    """
    ${base}
    """
    
}


process fasta_search
{
    errorStrategy 'ignore'
    
    input:
    val(id) from ch2
    
    output:
    val(id) into ch3
    
    script:
    def base_filename = id.toString().split("/")[-1]
    def base = """mkdir -p ${target_search_dir} && blastp -query ${target_fasta_dir}/${base_filename}.fa -db /root/Downloads/UP000005640_9606.fasta -out ${target_search_dir}/${base_filename}.blast-out.b6 -outfmt "6 std qseq sseq" -word_size 6 -evalue 0.05 -gapopen 11 -gapextend 1 -max_target_seqs 150"""
    """
    ${base}
    """
   
}

process filter_blast_results
{
    errorStrategy 'ignore'
    
    input:
    val(id) from ch3
    
    output:
    val(id) into ch4
    
    script:
    def base_filename = id.toString().split("/")[-1]
    def base = """mkdir -p ${target_search_dir_msa} && python ${filter_blast_script_path} ${target_fasta_dir}/${base_filename}.fa ${target_search_dir}/${base_filename}.blast-out.b6 ${homo_fasta_path} ${target_search_dir_msa}/${base_filename}.msa.fa"""
    """
    ${base}
    """
}


process align_with_mafft
{
    errorStrategy 'ignore'
    
    input:
    val(id) from ch4
    
    output:
    val(id) into ch5
    
    script:
    def base_filename = id.toString().split("/")[-1]
    def base = """mkdir -p ${target_search_dir_msa_out} && mafft ${target_search_dir_msa}/${base_filename}.msa.fa > ${target_search_dir_msa_out}/${base_filename}.msa.fa """
    """
    ${base}
    """
    
}
process final_step {
    input:
    val _ from ch5.collect()

    script:
    """
    cd ${target_search_dir_msa_out} && ${extract_script_path} > ${homology_count_output_path}
    """
}