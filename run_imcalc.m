function run_imcalc(input, output, directory, expression)

spm_jobman('initcfg');

matlabbatch{1}.spm.util.imcalc.input = input;
matlabbatch{1}.spm.util.imcalc.output = output;
matlabbatch{1}.spm.util.imcalc.outdir = {directory};
matlabbatch{1}.spm.util.imcalc.expression = expression; % if a mapis available, create a binarised map of non-zero values
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;

spm_jobman('run',matlabbatch);

end
