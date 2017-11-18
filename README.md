# ts-PFCTs
This is the project for functional feature extraction in fMRI registraion. We develop the tissue-specific patch-based functional correlation tensors (ts-PFCTs) to extract both gray-matter and white-matter functional information. 

Basic steps:

(1) input "mrData" and "maskData" corresponding to your preprocessed fMRI data and gray-matter/white matter mask

(2) output is named "tensor", this is the input of multi-channel registration.

For multi-channel registraion code, you can visit https://www.nitrc.org/projects/intergroupreg/ to download multi-channel Demons' code or email peizhang@email.unc.edu (Pei Zhang) to access multi-channel LDDMM's code

For future studies, you can also use FA, MD, RD, VR or other features maps as the input of multi-channel registration.

If you think it is useful for you, please star it. And you can cite this paper (this paper is a simplified version of our paper, a better version is under review by a journal):
"Zhou, Y., et al (2017). Improving Functional MRI Registration Using Whole-Brain Functional Correlation Tensors. MICCAI 2017, pp. 416-423."
