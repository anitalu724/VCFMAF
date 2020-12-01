vcf_list = ['examples/ms_vcf/colon1-sample.vcf',
'examples/ms_vcf/colon2-sample.vcf',
'examples/ms_vcf/colon3-sample.vcf',
'examples/ms_vcf/liver1-sample.vcf',
'examples/ms_vcf/liver2-sample.vcf',
'examples/ms_vcf/liver3-sample.vcf',
'examples/ms_vcf/intestine1-sample.vcf',
'examples/ms_vcf/intestine2-sample.vcf',
'examples/ms_vcf/intestine3-sample.vcf']

maf_list = [
'examples/ms_maf/colon1-sample.maf',
'examples/ms_maf/colon2-sample.maf',
'examples/ms_maf/colon3-sample.maf',
'examples/ms_maf/liver1-sample.maf',
'examples/ms_maf/liver2-sample.maf',
'examples/ms_maf/liver3-sample.maf',
'examples/ms_maf/intestine1-sample.maf',
'examples/ms_maf/intestine2-sample.maf',
'examples/ms_maf/intestine3-sample.maf']

vcf2vep2maf(vcf_list, maf_list, 'examples/ms_maf', )

