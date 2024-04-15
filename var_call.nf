#!/usr/bin/env nextflow

nextflow.enable.dsl=2

gene_name = params.gene
up_gene_symbol = gene_name.toUpperCase()
d_base = params.db_init
res_base = params.res_init
caller_base = params.caller_init
output_folder = params.out_dir

params.build='hg38'

if (params.build=='b37') {
    db = "${d_base}/${gene_name}/b37"
    res_dir = "${res_base}/${gene_name}/res_b37"
    caller_dir = "${caller_base}/${gene_name}/b37/bin"
    debug37 = "--minimum_extract_score_over_homref=0"
    debug38 = ""

    if (params.gene=='cyp2d6') {
        chrom = "22"
    	region_a1 = "22:42522000-42542000"
    	region_a2 = "042522000-042542000"
    	region_b1 = "22:42522300-42528400"
    	region_b2 = "042522300-042528400"

    } else if (params.gene=='cyp2a6') {
        chrom = "19"      
	region_a1 = "19:41339500-41389490"
	region_a2 = "041339500-041389490"
	region_b1 = "19:41348500-41358200"
	region_b2 = "041348500-041358200"

    } else if (params.gene=='cyp2b6') {	
        chrom = "19"
	region_a1 = "19:41487100-41534300"
	region_a2 = "041487100-041534300"
	region_b1 = "19:41494750-41522800"
	region_b2 = "041494750-041522800"

    } else if (params.gene=='cyp2c19') {
        chrom = "10"
	region_a1 = "10:96518000-96613000"
	region_a2 = "096518000-096613000"
	region_b1 = "10:96519000-96612750"
	region_b2 = "096519000-096612750"

    } else if (params.gene=='cyp2c9') {
        chrom = "10"
	region_a1 = "10:96695000-96749500"
	region_a2 = "096695000-096749500"
	region_b1 = "10:96695750-96748950"
	region_b2 = "096695750-096748950"

    } else if (params.gene=='cyp2c8') {
        chrom = "10"
	region_a1 = "10:96786500-96839250"
	region_a2 = "096786500-096839250"
	region_b1 = "10:96796750-96829600"
	region_b2 = "096796750-096829600"

    } else if (params.gene=='cyp2e1') {
        chrom = "10"
	region_a1 = "10:135330300-135352600"
	region_a2 = "135330300-135352600"
	region_b1 = "10:135339500-135351300"
	region_b2 = "135339500-135351300"

    } else if (params.gene=='cyp3a4') {
        chrom = "7"
	region_a1 = "7:99344600-99391800"
	region_a2 = "099344600-099391800"
	region_b1 = "7:99355450-99382600"
	region_b2 = "099355450-099382600"

    } else if (params.gene=='cyp3a5') {
        chrom = "7"
	region_a1 = "7:99235800-99287600"
	region_a2 = "099235800-099287600"
	region_b1 = "7:99245850-99278250"
	region_b2 = "099245850-099278250"

    } else if (params.gene=='cyp1a1') {
        chrom = "15"
	region_a1 = "15:75001900-75027900"
	region_a2 = "075001900-075027900"
	region_b1 = "15:75011600-75015320"
	region_b2 = "075011600-075015320"

    } else if (params.gene=='cyp1a2') {
        chrom = "15"
        region_a1 = "15:75031200-75058900"
        region_a2 = "075031200-075058900"
        region_b1 = "15:75038200-75047650"
        region_b2 = "075038200-075047650"

    } else if (params.gene=='cyp4f2') {
        chrom = "19"
	region_a1 = "19:15978850-16018850"
	region_a2 = "015978850-016018850"
	region_b1 = "19:15988850-16008850"
	region_b2 = "015988850-016008850"

    } else if (params.gene=='cypor') {
       chrom = "7"
       region_a1 = "7:75534473-75626173"
       region_a2 = "075534473-075626173"
       region_b1 = "7:75549473-75621173"
       region_b2 = "075549473-075621173"
    }


} else if (params.build=='hg19') {
    db = "${d_base}/${gene_name}/b37"
    res_dir = "${res_base}/${gene_name}/res_hg19"
    caller_dir = "${caller_base}/${gene_name}/b37/bin"
    debug37 = "--minimum_extract_score_over_homref=0"
    debug38 = ""

    if (params.gene=='cyp2d6') {
        chrom = "chr22"
    	region_a1 = "chr22:42522000-42542000"
    	region_a2 = "042522000-042542000"
    	region_b1 = "chr22:42522300-42528400"
    	region_b2 = "042522300-042528400"

    } else if (params.gene=='cyp2a6') {
        chrom = "chr19"
        region_a1 = "chr19:41339500-41389490"
        region_a2 = "041339500-041389490"
        region_b1 = "chr19:41348500-41358200"
        region_b2 = "041348500-041358200"


    } else if (params.gene=='cyp2b6') {
      	chrom = "chr19"
        region_a1 = "chr19:41487100-41534300"
        region_a2 = "041487100-041534300"
	region_b1 = "chr19:41494750-41522800"
	region_b2 = "041494750-041522800"

    } else if (params.gene=='cyp2c19') {
        chrom = "chr10"
         region_a1 = "chr10:96518000-96613000"
         region_a2 = "096518000-096613000"
         region_b1 = "chr10:96519000-96612750"
         region_b2 = "096519000-096612750"

    } else if (params.gene=='cyp2c9') {
        chrom = "chr10"
        region_a1 = "chr10:96695000-96749500"
        region_a2 = "096695000-096749500"
        region_b1 = "chr10:96695750-96748950"
        region_b2 = "096695750-096748950"

    } else if (params.gene=='cyp2c8') {
        chrom = "chr10"
	region_a1 = "chr10:96786500-96839250"
	region_a2 = "096786500-096839250"
	region_b1 = "chr10:96796750-96829600"
	region_b2 = "096796750-096829600"


    } else if (params.gene=='cyp2e1') {
        chrom = "chr10"
        region_a1 = "chr10:135330300-135352600"
	region_a2 = "135330300-135352600"
	region_b1 = "chr10:135339500-135351300"
	region_b2 = "135339500-135351300"

    } else if (params.gene=='cyp3a4') {
      	chrom = "chr7"
	region_a1 = "chr7:99344600-99391800"
	region_a2 = "099344600-099391800"
	region_b1 = "chr7:99355450-99382600"
	region_b2 = "099355450-099382600"


    } else if (params.gene=='cyp3a5') {
      	chrom = "chr7"
	region_a1 = "chr7:99235800-99287600"
	region_a2 = "099235800-099287600"
	region_b1 = "chr7:99245850-99278250"
	region_b2 = "099245850-099278250"

    } else if (params.gene=='cyp1a1') {
      	chrom = "chr15"
        region_a1 = "chr15:75001900-75027900"
	region_a2 = "075001900-075027900"
	region_b1 = "chr15:75011600-75015320"
	region_b2 = "075011600-075015320"


    } else if (params.gene=='cyp1a2') {
        chrom = "chr15"
	region_a1 = "chr15:75031200-75058900"
    	region_a2 = "075031200-075058900"
	region_b1 = "chr15:75038200-75047650"
	region_b2 = "075038200-075047650"

    } else if (params.gene=='cyp4f2') {
        chrom = "chr19"
        region_a1 = "chr19:15978850-16018850"
	region_a2 = "015978850-016018850"
	region_b1 = "chr19:15988850-16008850"
	region_b2 = "015988850-016008850"

    } else if (params.gene=='cypor') {
       chrom = "chr7"
       region_a1 = "chr7:75534473-75626173"
       region_a2 = "075534473-075626173"
       region_b1 = "chr7:75549473-75621173"
       region_b2 = "075549473-075621173"

    }


} else {
    db = "${d_base}/${gene_name}/hg38"
    res_dir = "${res_base}/${gene_name}/res_hg38"
    caller_dir = "${caller_base}/${gene_name}/hg38/bin"
    debug37 = ""
    debug38 = "--minimum_extract_score_over_homref=0"

    if (params.gene=='cyp2d6') {
        chrom = "chr22"
    	region_a1 = "chr22:42126000-42137500"
    	region_a2 = "042126000-042137500"
    	region_b1 = "chr22:42126300-42132400"
    	region_b2 = "042126300-042132400"
	transcript = "ENST00000645361"

    } else if (params.gene=='cyp2a6') {
        chrom = "chr19"
	region_a1 = "chr19:40833541-40887450"
	region_a2 = "040833541-040887450"
	region_b1 = "chr19:40842750-40852250"
	region_b2 = "040842750-040852250"
	transcript = "ENST00000301141"
    
    } else if (params.gene=='cyp2b6') {
        chrom = "chr19"
	region_a1 = "chr19:40981280-41028400"
	region_a2 = "040981280-041028400"
	region_b1 = "chr19:40988800-41016900"
	region_b2 = "040988800-041016900"
	transcript = "ENST00000324071"

    } else if (params.gene=='cyp2c19') {
        chrom = "chr10"
	region_a1 = "chr10:94752750-94865500"
	region_a2 = "094752750-094865500"
	region_b1 = "chr10:94759250-94853000"
	region_b2 = "094759250-094853000"
	transcript = "ENST00000371321"

    } else if (params.gene=='cyp2c9') {
      	chrom = "chr10"
	region_a1 = "chr10:94935000-94990000"
	region_a2 = "094935000-094990000"
	region_b1 = "chr10:94936000-94989200"
	region_b2 = "094936000-094989200"
	transcript = "ENST00000260682"

    } else if (params.gene=='cyp2c8') {
        chrom = "chr10"
	region_a1 = "chr10:95026750-95079500"
	region_a2 = "095026750-095079500"
	region_b1 = "chr10:95037050-95069800"
	region_b2 = "095037050-095069800"	
	transcript = "ENST00000371270"

    } else if (params.gene=='cyp2e1') {
        chrom = "chr10"
	region_a1 = "chr10:133517350-133549100"
	region_a2 = "133517350-133549100"
	region_b1 = "chr10:133522363-133544123"
	region_b2 = "133522363-133544123"
	transcript = "ENST00000252945"

    } else if (params.gene=='cyp3a4') {
        chrom = "chr7"
	region_a1 = "chr7:99746950-99794100"
	region_a2 = "099746950-099794100"
	region_b1 = "chr7:99757850-99784950"
	region_b2 = "099757850-099784950"
	transcript = "ENST00000651514"

    } else if (params.gene=='cyp3a5') {	
        chrom = "chr7"
	region_a1 = "chr7:99638200-99690000"
	region_a2 = "099638200-099690000"
	region_b1 = "chr7:99648250-99680650"	
	region_b2 = "099648250-099680650"
	transcript = "ENST00000222982"

    } else if (params.gene=='cyp1a1') {
        chrom = "chr15"
	region_a1 = "chr15:74709550-74735500"
	region_a2 = "074709550-074735500"
	region_b1 = "chr15:74719250-74725500"
	region_b2 = "074719250-074725500"
	transcript = "ENST00000379727"

    } else if (params.gene=='cyp1a2') {
        chrom = "chr15"
	region_a1 = "chr15:74738850-74766600"
	region_a2 = "074738850-074766600"
	region_b1 = "chr15:74745850-74755280"
	region_b2 = "074745850-074755280"
	transcript = "ENST00000343932"

    } else if (params.gene=='cyp4f2') {
        chrom = "chr19"
	region_a1 = "chr19:15868100-15908000"
	region_a2 = "015868100-015908000"
	region_b1 = "chr19:15878000-15898100"
	region_b2 = "015878000-015898100"
	transcript = "ENST00000221700"

    } else if (params.gene=='cypor') {
       chrom = "chr7"
       region_a1 = "chr7:75905155-75996855"
       region_a2 = "075905155-075996855"
       region_b1 = "chr7:75910155-75991855"
       region_b2 = "075910155-075991855"
       transcript = "ENST00000461988"

    } else if (params.gene=='nat1') {
       chrom = "chr8"
       region_a1 = "chr8:18170000-18226000"
       region_a2 = "018170000-018226000"
       region_b1 = "chr8:18221000-18224000"
       region_b2 = "018221000-018224000"
       transcript = "ENST00000307719"

    } else if (params.gene=='nat2') {
       chrom = "chr8"
       region_a1 = "chr8:18386000-18402000"
       region_a2 = "018386000-018402000"
       region_b1 = "chr8:18386000-18402000"
       region_b2 = "018386000-018402000"
       transcript = "ENST00000286479"

    } else if (params.gene=='gstm1') {
       chrom = "chr1"
       region_a1 = "chr1:109687000-109695000"
       region_a2 = "109687000-109695000"
       region_b1 = "chr1:109687000-109695000"
       region_b2 = "109687000-109695000"
       transcript = "ENST00000309851"

    } else if (params.gene=='gstt1') {
       chrom = "chr22_KI270879v1_alt"
       region_a1 = "chr22_KI270879v1_alt:269000-280000"
       region_a2 = "000269000-000280000"
       region_b1 = "chr22_KI270879v1_alt:269000-280000"
       region_b2 = "000269000-000280000"
       transcript = "ENST00000612885"

    }  else if (params.gene=='slco1b1') {
       chrom = "chr12"
       region_a1 = "chr12:21120000-21250000"
       region_a2 = "021120000-021250000"
       region_b1 = "chr12:21120000-21250000"
       region_b2 = "021120000-021250000"
       transcript = "ENST00000256958"
       
    }  else if (params.gene=='nudt15') {
       chrom = "chr13"
       region_a1 = "chr13:48037500-48055000"
       region_a2 = "048037500-048055000"
       region_b1 = "chr13:48037500-48055000"
       region_b2 = "048037500-048055000"
       transcript = "ENST00000258662"

    }  else if (params.gene=='tpmt') {
       chrom = "chr6"
       region_a1 = "chr6:18128000-18155500"
       region_a2 = "018128000-018155500"
       region_b1 = "chr6:18128500-18155000"
       region_b2 = "018128500-018155000"
       transcript = "ENST00000309983"

    }  else if (params.gene=='ugt1a1') {
       chrom = "chr2"
       region_a1 = "chr2:233758500-233772600"
       region_a2 = "233758500-233772600"
       region_b1 = "chr2:233759000-233772000"
       region_b2 = "233759000-233772000"
       transcript = "ENST00000305208"

    }  else if (params.gene=='dpyd') {
       chrom = "chr1"
       region_a1 = "chr1:97072500-97926000"
       region_a2 = "097072500-097926000"
       region_b1 = "chr1:97072500-97926000"
       region_b2 = "097072500-097926000"
       transcript = "ENST00000370192"
    }
}

params.format='binary'

if (params.format=='compressed') {
    ext = "cram"
    ind = "crai"
    cram_options = "--force_use_input_ref_for_cram_reading" 

} else { 
    ext = "bam"
    ind = "bai"
    cram_options = ""

}


ref_dir_val = new File("${params.ref_file}").getParent()
ref_genome = new File("${params.ref_file}").getName()

workflow {
    ref_ch = Channel.value("${ref_dir_val}")
    data_ch = Channel.fromFilePairs(params.in_bam, type: 'file') {  file -> file.name.replaceAll(/.${ext}|.${ind}$/,'') }
    call_snvs1(data_ch, ref_ch, res_dir)
    call_snvs2(data_ch, ref_ch)
    get_depth(data_ch, ref_ch, res_dir)
    call_snvs1.out.join(call_snvs2.out).set {var_ch_joined}
    format_snvs(var_ch_joined)
    save_vcfs(format_snvs.out)
    get_core_var(format_snvs.out, res_dir, res_base, ref_ch, caller_base)
    get_core_var.out.join(get_depth.out).set {fin_files}
    analyse(fin_files, caller_dir)
}


process call_snvs1 {
//   maxForks 10

    input:
    tuple val(name), path(bam)
    path ref_dir
    path res_dir

    output:	        
    tuple val(name), path("${name}_var_1")
      
    script:
    ext1 = bam[0].getExtension()
    if (ext1=='bam')
        sam_ind='bam.bai'

    else if (ext1=='cram')
        sam_ind='cram.crai'

    else if (ext1=='crai')
        sam_ind='crai'

    else
        sam_ind='bai'   
      
    """
    graphtyper genotype ${ref_dir}/${ref_genome} --sam=${name}.${ext} --sams_index=<(echo ${name}.${sam_ind}) --region=${region_a1} --output=${name}_var_1 --prior_vcf=${res_dir}/common_plus_core_var.vcf.gz -a ${debug38} ${cram_options}
    bcftools concat ${name}_var_1/${chrom}/*.vcf.gz > ${name}_var_1/${chrom}/${region_a2}.vcf 
    bgzip -f ${name}_var_1/${chrom}/${region_a2}.vcf 
    tabix -f ${name}_var_1/${chrom}/${region_a2}.vcf.gz

    """

}


process call_snvs2 {
//   maxForks 10

    input:
    tuple val(name), path(bam)
    path ref_dir

    output: 
    tuple val(name), path("${name}_var_2")

    script:
    ext1 = bam[0].getExtension()
    if (ext1=='bam')
        sam_ind='bam.bai'

    else if (ext1=='cram')
        sam_ind='cram.crai'

    else if (ext1=='crai')
        sam_ind='crai'

    else 
        sam_ind='bai'

    """
    graphtyper genotype ${ref_dir}/${ref_genome} --sam=${name}.${ext} --sams_index=<(echo ${name}.${sam_ind}) --region=${region_a1} --output=${name}_var_2 -a ${debug38} ${debug37} ${cram_options}
    bcftools concat ${name}_var_2/${chrom}/*.vcf.gz > ${name}_var_2/${chrom}/${region_a2}.vcf       
    bgzip -f ${name}_var_2/${chrom}/${region_a2}.vcf
    tabix -f ${name}_var_2/${chrom}/${region_a2}.vcf.gz

    """

}


process get_depth {
//   maxForks 10

    input:
    tuple val(name), path(bam)
    path ref_dir
    path res_dir

    output:
    tuple val(name), path("${name}_${gene_name}_ctrl.depth")

    script:

    """   
    samtools bedcov --reference ${ref_dir}/${ref_genome} ${res_dir}/test3.bed ${name}.${ext} > ${name}_${gene_name}_ctrl.depth      

    """

}


process format_snvs {
//   maxForks 10

//    publishDir "$output_folder/$gene_name/variants", mode: 'copy', overwrite: 'true', saveAs: { ${name}_var/${name}_all_norm.vcf.gz -> "${name}_${gene_name}.vcf.gz" }

    input:
    tuple val(name), path("${name}_var_1"), path("${name}_var_2")

    output:
    tuple val(name), path("*")

    script:

    """
    bcftools isec -p ${name}_var -Oz ${name}_var_1/${chrom}/${region_a2}.vcf.gz ${name}_var_2/${chrom}/${region_a2}.vcf.gz
    bcftools concat -a -D -r ${region_b1} ${name}_var/0000.vcf.gz ${name}_var/0001.vcf.gz ${name}_var/0002.vcf.gz -Oz -o ${name}_var/${name}_${region_b2}.vcf.gz
    tabix ${name}_var/${name}_${region_b2}.vcf.gz
    bcftools norm -m - ${name}_var/${name}_${region_b2}.vcf.gz | bcftools view -e 'GT="1/0"' | bcftools view -e 'GT="0/0"' | bcftools view -e 'FILTER="PASS" & INFO/QD<10 || 0<ABHet<0.25' | bgzip -c > ${name}_var/${name}_all_norm.vcf.gz
    tabix ${name}_var/${name}_all_norm.vcf.gz

    """

}


process save_vcfs {

    publishDir "$output_folder/$gene_name/variants", pattern: '*vcf.gz', mode: 'copy', overwrite: 'true'
    publishDir "$output_folder/$gene_name/variants", pattern: '*vcf.gz.tbi', mode: 'copy', overwrite: 'true'

    errorStrategy 'ignore'
    tag "${name}"

    input:
    tuple val(name), path("${name}_vars")

    output:
    tuple val(name), path("*")

    script:
    """
    cp ${name}_vars/${name}_all_norm.vcf.gz ./${name}_${gene_name}.vcf.gz
    cp ${name}_vars/${name}_all_norm.vcf.gz.tbi ./${name}_${gene_name}.vcf.gz.tbi
    """

}

process get_core_var {
//   maxForks 10
   
    errorStrategy 'ignore'
    tag "${name}"   

    input:
    tuple val(name), path("${name}_vars")
    path res_dir
    path res_base
    path ref_dir
    path caller_base

    output:
    tuple val(name), path("${name}_int")

    script:
 
    if (params.build=='b37') {

    """
    bcftools isec ${name}_vars/${name}_all_norm.vcf.gz ${res_dir}/allele_def_var.vcf.gz -p ${name}_int -Oz
    bcftools norm -m - ${name}_int/0002.vcf.gz | bcftools view -e 'GT="1/0"' | bcftools view -e 'GT="0/0"' | bgzip -c > ${name}_int/${name}_core.vcf.gz
    tabix ${name}_int/${name}_core.vcf.gz
    """

    } else if (params.build=='hg19') {

    """
    bcftools isec ${name}_vars/${name}_all_norm.vcf.gz ${res_dir}/allele_def_var.vcf.gz -p ${name}_int -Oz
    bcftools norm -m - ${name}_int/0002.vcf.gz | bcftools view -e 'GT="1/0"' | bcftools view -e 'GT="0/0"' | bgzip -c > ${name}_int/${name}_core.vcf.gz
    tabix ${name}_int/${name}_core.vcf.gz
    """

    } else {

    """
    bcftools csq -p m -v 0 -f ${ref_dir}/${ref_genome} -g ${res_base}/annotation/Homo_sapiens.GRCh38.110.gff3.gz ${name}_vars/${name}_all_norm.vcf.gz -o ${name}_vars/${name}_all_norm_annot.vcf
    bgzip ${name}_vars/${name}_all_norm_annot.vcf
    tabix ${name}_vars/${name}_all_norm_annot.vcf.gz

    bcftools isec ${name}_vars/${name}_all_norm_annot.vcf.gz ${res_dir}/allele_def_var.vcf.gz -p ${name}_int -Oz
    bcftools norm -m - ${name}_int/0002.vcf.gz | bcftools view -e 'GT="1/0"' | bcftools view -e 'GT="0/0"' > ${name}_int/${name}_core_int1.vcf
    
    bgzip -d ${name}_int/0000.vcf.gz
    python3 ${caller_base}/novel/core_var.py ${name}_int/0000.vcf ${up_gene_symbol} ${transcript} >> ${name}_int/${name}_core_int1.vcf
    bcftools sort ${name}_int/${name}_core_int1.vcf -T ${name}_int | bgzip -c > ${name}_int/${name}_core.vcf.gz
    tabix ${name}_int/${name}_core.vcf.gz

    """
    }
}



process analyse {

    publishDir "$output_folder/$gene_name/output", mode: 'copy', overwrite: 'true'

    errorStrategy 'ignore'
    tag "${name}"

    input:
    tuple val(name), path("${name}_int"), path("${name}_${gene_name}_dp")
    path caller_dir

    output:
    tuple val(name), file("${name}_${gene_name}.info")

    script:
   
    """
    bcftools query -f'[%POS~%REF>%ALT~%GT\n]' ${name}_int/${name}_core.vcf.gz > ${name}_core_snvs.dip
    python3 ${caller_dir}/stellarpgx.py ${name}_core_snvs.dip ${name}_${gene_name}_dp > ${name}_${gene_name}.info

    """

}

