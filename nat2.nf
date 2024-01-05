#!/usr/bin/env nextflow

nextflow.enable.dsl=1

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

    } else if (params.gene=='cyp2a6') {
        chrom = "chr19"
	region_a1 = "chr19:40833541-40887450"
	region_a2 = "040833541-040887450"
	region_b1 = "chr19:40842750-40852250"
	region_b2 = "040842750-040852250"
    
    } else if (params.gene=='cyp2b6') {
        chrom = "chr19"
	region_a1 = "chr19:40981280-41028400"
	region_a2 = "040981280-041028400"
	region_b1 = "chr19:40988800-41016900"
	region_b2 = "040988800-041016900"

    } else if (params.gene=='cyp2c19') {
        chrom = "chr10"
	region_a1 = "chr10:94752750-94865500"
	region_a2 = "094752750-094865500"
	region_b1 = "chr10:94759250-94853000"
	region_b2 = "094759250-094853000"

    } else if (params.gene=='cyp2c9') {
      	chrom = "chr10"
	region_a1 = "chr10:94935000-94990000"
	region_a2 = "094935000-094990000"
	region_b1 = "chr10:94936000-94989200"
	region_b2 = "094936000-094989200"

    } else if (params.gene=='cyp2c8') {
        chrom = "chr10"
	region_a1 = "chr10:95026750-95079500"
	region_a2 = "095026750-095079500"
	region_b1 = "chr10:95037050-95069800"
	region_b2 = "095037050-095069800"	

    } else if (params.gene=='cyp2e1') {
        chrom = "chr10"
	region_a1 = "chr10:133517350-133549100"
	region_a2 = "133517350-133549100"
	region_b1 = "chr10:133526050-133537800"
	region_b2 = "133526050-133537800"

    } else if (params.gene=='cyp3a4') {
        chrom = "chr7"
	region_a1 = "chr7:99746950-99794100"
	region_a2 = "099746950-099794100"
	region_b1 = "chr7:99757850-99784950"
	region_b2 = "099757850-099784950"

    } else if (params.gene=='cyp3a5') {	
        chrom = "chr7"
	region_a1 = "chr7:99638200-99690000"
	region_a2 = "099638200-099690000"
	region_b1 = "chr7:99648250-99680650"	
	region_b2 = "099648250-099680650"

    } else if (params.gene=='cyp1a1') {
        chrom = "chr15"
	region_a1 = "chr15:74709550-74735500"
	region_a2 = "074709550-074735500"
	region_b1 = "chr15:74719250-74725500"
	region_b2 = "074719250-074725500"

    } else if (params.gene=='cyp1a2') {
        chrom = "chr15"
	region_a1 = "chr15:74738850-74766600"
	region_a2 = "074738850-074766600"
	region_b1 = "chr15:74745850-74755280"
	region_b2 = "074745850-074755280"


    } else if (params.gene=='cyp4f2') {
        chrom = "chr19"
	region_a1 = "chr19:15868100-15908000"
	region_a2 = "015868100-015908000"
	region_b1 = "chr19:15878000-15898100"
	region_b2 = "015878000-015898100"

    } else if (params.gene=='cypor') {
       chrom = "chr7"
       region_a1 = "chr7:75905155-75996855"
       region_a2 = "075905155-075996855"
       region_b1 = "chr7:75910155-75991855"
       region_b2 = "075910155-075991855"

    } else if (params.gene=='nat2') {
       chrom = "chr8"
       region_a1 = "chr8:18386000-18402000"
       region_a2 = "018386000-018402000"
       region_b1 = "chr8:18386000-18402000"
       region_b2 = "018386000-018402000"

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


align_file = Channel.fromFilePairs(params.in_bam, type: 'file') {  file -> file.name.replaceAll(/.${ext}|.${ind}$/,'') }

align_file.into { data1; data2; data3; data4; data5; data6 }


ref_dir_val = new File("${params.ref_file}").getParent()
ref_genome = new File("${params.ref_file}").getName()



process call_snvs1 {
//   maxForks 10

    input:
    set val(name), file(bam) from data1
    path ref_dir from Channel.value("${ref_dir_val}")
    path res_dir

    output:	        
    set val(name), path("${name}_var_1") into var_ch1
      
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
    set val(name), file(bam) from data2
    path ref_dir from Channel.value("${ref_dir_val}")

    output: 
    set val(name), path("${name}_var_2") into var_ch2

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


process call_sv_del {
//   maxForks 10

    input:
    set val(name), file(bam) from data3
    path ref_dir from Channel.value("${ref_dir_val}")
    path res_dir

    output:
    set val(name), path("${name}_sv_del") into sv_ch1

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
    graphtyper genotype_sv ${ref_dir}/${ref_genome} --sam=${name}.${ext} --region=${region_a1} --output=${name}_sv_del ${res_dir}/sv_test.vcf.gz

    """

}

process call_sv_dup {
//   maxForks 10

    input:
    set val(name), file(bam) from data4
    path ref_dir from Channel.value("${ref_dir_val}")
    path res_dir

    output:
    set val(name), path("${name}_sv_dup") into sv_ch2

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
    graphtyper genotype_sv ${ref_dir}/${ref_genome} --sam=${name}.${ext} --region=${region_a1} --output=${name}_sv_dup ${res_dir}/sv_test3.vcf.gz

    """
}


process get_depth {
//   maxForks 10

    input:
    set val(name), file(bam) from data5
    path ref_dir from Channel.value("${ref_dir_val}")
    path res_dir

    output:
    set val(name), file("${name}_${gene_name}_ctrl.depth") into sv_ch3

    script:

    """   
    samtools bedcov --reference ${ref_dir}/${ref_genome} ${res_dir}/test3.bed ${name}.${ext} > ${name}_${gene_name}_ctrl.depth      

    """

}


var_ch1.join(var_ch2).set { var_ch_joined }


process format_snvs {
//   maxForks 10

    publishDir "$output_folder/$gene_name/variants", pattern: '*vcf.gz', mode: 'copy', overwrite: 'true'
    publishDir "$output_folder/$gene_name/variants", pattern: '*vcf.gz.tbi', mode: 'copy', overwrite: 'true'    

    input:
    set val(name), path("${name}_var_1"), path("${name}_var_2") from var_ch_joined

    output:
    set val(name), path("${name}_var") into (var_norm1, var_norm2)
    set val(name), file("${name}_${gene_name}.vcf.gz"), file("${name}_${gene_name}.vcf.gz.tbi") into var_norm3

    script:

    """
        bcftools isec -p ${name}_var -Oz ${name}_var_1/${chrom}/${region_a2}.vcf.gz ${name}_var_2/${chrom}/${region_a2}.vcf.gz
        bcftools concat -a -D -r ${region_b1} ${name}_var/0000.vcf.gz ${name}_var/0001.vcf.gz ${name}_var/0002.vcf.gz -Oz -o ${name}_var/${name}_${region_b2}.vcf.gz
        tabix ${name}_var/${name}_${region_b2}.vcf.gz
        bcftools norm -m - ${name}_var/${name}_${region_b2}.vcf.gz | bcftools view -e 'GT="1/0"' | bcftools view -e 'GT="0/0"' | bcftools view -e 'FILTER="PASS" & INFO/QD<10 || 0<ABHet<0.25' | bgzip -c > ${name}_var/${name}_all_norm.vcf.gz
        tabix ${name}_var/${name}_all_norm.vcf.gz
	cp ${name}_var/${name}_all_norm.vcf.gz ./${name}_${gene_name}.vcf.gz
	cp ${name}_var/${name}_all_norm.vcf.gz.tbi ./${name}_${gene_name}.vcf.gz.tbi

    """

}


process get_core_var {
//   maxForks 10
   
    errorStrategy 'ignore'
    tag "${name}"   

    input:
    set val(name), path("${name}_vars") from var_norm1
    path res_dir
    path ref_dir from Channel.value("${ref_dir_val}")

    output:
    set val(name), path("${name}_int") into (core_vars1, core_vars2, core_vars3)

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
    bcftools isec ${name}_vars/${name}_all_norm.vcf.gz ${res_dir}/allele_def_var.vcf.gz -p ${name}_int -Oz
    bcftools norm -m - ${name}_int/0002.vcf.gz | bcftools view -e 'GT="1/0"' | bcftools view -e 'GT="0/0"' > ${name}_int/${name}_core_int1.vcf
    bcftools csq -p m -v 0 -f ${ref_dir}/${ref_genome} -g ${res_base}/annotation/Homo_sapiens.GRCh38.110.gff3.gz ${name}_int/0000.vcf.gz -o ${name}_int/0000_annot.vcf

    python3 ${caller_base}/novel/core_var.py ${name}_int/0000_annot.vcf ${up_gene_symbol} ${transcript} >> ${name}_int/${name}_core_int1.vcf
    bcftools sort ${name}_int/${name}_core_int1.vcf -T ${name}_int | bgzip -c > ${name}_int/${name}_core.vcf.gz
    tabix ${name}_int/${name}_core.vcf.gz

    """
    }

}


core_vars3.join(data6)
.set { phase_ch }
    
process phase_snvs {
    label 'phase'    

    errorStrategy 'ignore'
    tag "${name}"

    // publishDir "$output_folder/phased_vcfs", mode: 'copy', overwrite: 'true'

    input:
    set val(name), path("${name}_int1"), file(bam) from phase_ch
    path ref_dir from Channel.value("${ref_dir_val}")

    output:
    set val(name), file("${name}_${gene_name}_phased.vcf.gz") into phased_ch

    script:

    """
    whatshap phase --output ${name}_${gene_name}_phased.vcf.gz --reference  ${ref_dir}/${ref_genome} --chromosome $chrom ${name}_int1/${name}_core.vcf.gz ${name}.${ext} --indels
        
    """
}


process analyse_1 {
//   maxForks 10

    errorStrategy 'ignore'
    tag "${name}"

    input:
    set val(name), path("${name}_gene_del") from sv_ch1

    output:
    set val(name), path("${name}_gene_del/${name}_gene_del_summary.txt") into del_ch

    script:

    """
    bcftools query -f'%ID\t%ALT\t[%GT\t%DP]\t%INFO/ABHet\t%INFO/ABHom\n' ${name}_gene_del/${chrom}/${region_a2}.vcf.gz > ${name}_gene_del/${name}_gene_del_summary.txt
    """

}


sv_ch2.join(core_vars1).set {dup_int}

process analyse_2 {
//   maxForks 10

    errorStrategy 'ignore'
    tag "${name}"

    input:
    set val(name), path("${name}_gene_dup"), path("${name}_int") from dup_int

    output:
    set val(name), path("${name}_gene_dup/${name}_gene_dup_summary.txt") into dup_ch

    script:

    """
    bcftools query -f'%POS~%REF>%ALT\t[%GT\t%DP]\t%INFO/ABHet\t%INFO/ABHom\n' -i'GT="alt"' ${name}_gene_dup/${chrom}/${region_a2}.vcf.gz > ${name}_gene_dup/${name}_gene_dup_summary.txt
    bcftools query -f'%POS~%REF>%ALT\t[%GT\t%DP]\t%INFO/ABHet\t%INFO/ABHom\n' -i'GT="alt"' ${name}_int/${name}_core.vcf.gz >> ${name}_gene_dup/${name}_gene_dup_summary.txt

    """

}


var_norm2.join(core_vars2).set {dip_req}

process analyse_3 {
//   maxForks 10

    errorStrategy 'ignore'
    tag "${name}"

    input:
    set val(name), path("${name}_vars"), path("${name}_int") from dip_req

    output:
    set val(name), path("${name}_vars/${name}_core_snvs.dip"), path("${name}_vars/${name}_full.dip"), path("${name}_vars/${name}_gt.dip") into prep_ch

    script:
    """
    bcftools query -f'[%POS~%REF>%ALT~%GT\n]' ${name}_int/${name}_core.vcf.gz  > ${name}_vars/${name}_core_snvs.dip
    bcftools query -f '%POS~%REF>%ALT\n' ${name}_vars/${name}_all_norm.vcf.gz > ${name}_vars/${name}_full.dip
    bcftools query -f'[%POS~%REF>%ALT~%GT\n]' ${name}_vars/${name}_all_norm.vcf.gz > ${name}_vars/${name}_gt.dip

    """

}


process analyse_4 {

    errorStrategy 'ignore'
    tag "${name}"

    input:
    set val(name), file(phased_vcf) from phased_ch

    output:
    set val(name), file("${name}_${gene_name}_hap_snvs.dip") into hap_snvs_ch

    script:
    """
    tabix ${phased_vcf}
    bcftools query -e 'GT="0|1" || GT="0/1"' -f'[%POS~%REF>%ALT;]' ${phased_vcf} | rev | cut -d";" -f2- | rev > ${name}_${gene_name}_hap_snvs.dip
    bcftools query -e 'GT="1|0" || GT="0/1"' -f'[%POS~%REF>%ALT;]' ${phased_vcf} | rev | cut -d";" -f2- | rev >> ${name}_${gene_name}_hap_snvs.dip
    """


}


prep_ch.join(del_ch).set {fin_files1}
fin_files1.join(dup_ch).set {fin_files2}
fin_files2.join(sv_ch3).set {fin_files3}
fin_files3.join(hap_snvs_ch).set {fin_files}


process call_stars {
//   maxForks 10

    publishDir "$output_folder/$gene_name/alleles", mode: 'copy', overwrite: 'true'

    errorStrategy 'ignore'
    tag "${name}"

    input:
    set val(name), path("${name}_vars/${name}_core_snvs.dip"), path("${name}_vars/${name}_full.dip"), path("${name}_vars/${name}_gt.dip"), path("${name}_gene_del/${name}_gene_del_summary.txt"), path("${name}_gene_dup/${name}_gene_dup_summary.txt"), file("${name}_${gene_name}_dp"), file("${name}_${gene_name}_hap_snvs.dip") from fin_files
    path db
    path caller_dir

    output:
    set val(name), file("${name}_${gene_name}.alleles") into star_ch

    script:
   
    """
    python3 ${caller_dir}/stellarpgx.py ${db}/diplo_db_debugged2.dbs ${name}_vars/${name}_core_snvs.dip ${name}_vars/${name}_full.dip ${name}_vars/${name}_gt.dip ${db}/genotypes4.dbs ${name}_gene_del/${name}_gene_del_summary.txt ${name}_gene_dup/${name}_gene_dup_summary.txt ${name}_${gene_name}_dp ${name}_${gene_name}_hap_snvs.dip ${db}/haps_var_new.dbs ${db}/a_scores.dbs > ${name}_${gene_name}.alleles  

    """

}

