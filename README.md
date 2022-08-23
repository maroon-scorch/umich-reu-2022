# umich-reu-2022
This repository hosts some of the miscellaneous codes and writeups I wrote for the [University of Michigan's 2022 Mathematics REU program](https://lsa.umich.edu/math/undergraduates/research-and-career-opportunities/research/research-experience-for-undergraduates--reu-.html) advised by Dr. [Lena Ji](http://www-personal.umich.edu/~lenaji/).

If you would like to learn more about our project, feel free to check out the paper [''Rationality of real conic bundles with quartic discriminant curve''](https://arxiv.org/abs/2208.08916) by Lena Ji and Mattie Ji. For a softer introduction to the subject, check out my presentation for the [UMich REU Seminar](Presentation/UMich_REU_Presentation/UMich_Reu_Presentation.pdf) and at the [Young Mathematician's Conference](Presentation/YMC_Presentation_Summer_2022/YMC_Presentation_Summer_2022.pdf).

My REU was funded by the National Science Foundation (Karen Smith’s NSF grant DMS-2101075).

## Important Note

1. Some of my **Sage** code relied on the **Sage** [code](http://sites.math.washington.edu/~vinzant/research/quartics/quartictype.sage) written by Daniel Plaumann, Bernd Sturmfels, and Cynthia Vinzant as supplementary material for their paper ["Quartic Curves and Their Bitangents"](https://arxiv.org/abs/1008.4104). Specifically, any functionality of the following list used their code
    - (a) checked for the topological type of a given smooth quartic
    - (b) found the bitangents of a given smooth quartic
    - (c) checked if $\Delta$ is smooth, used their code.
Any example listed with a known topological type was confirmed by their code as well.

I acknowledge this fact. Therefore, I have replaced ALL the places where I used the code in Case (a) and (b) with blank functions instead. If you are curious on what the code should have looked like with their code, feel free to check out the accomapnied [repository](https://github.com/lena-ji/ConicBundles) for our paper or try to insert Plaumann et al.'s original code.

The only exception is in Case (c), this is because Plaumann et al.'s smoothness check is a 2-line code check using the Jacobian Criterion for Smoothness, which is a well-known fact in Algebraic Geometry (see pg. 31 of Hartshorne's **Algebraic Geometry**). Since this is general enough, I have kept it in.

2. Some of my **Sage** code were based on the [Magma code](https://github.com/ivogt161/FJSVV-rationality) accompanying the paper ["Curve classes on conic bundle threefolds and applications to rationality"](https://arxiv.org/abs/2207.07093) by Sarah Frei, Lena Ji, Soumya Sankar, Bianca Viray, and Isabel Vogt. The part in **Sage** code that referenced their work were re-implementations of functionalities from their **Magma** code. Specifically, any part where I checked the signatures of the fibers of $\pi_1$ were taken from or inspired by their code. The actual file where I wrote the re-implemenation is [here](Code/connected%20components/general_qts.sage), any other **Sage** file that references the work of Frei et al. will just use the code from the aforementioned file.

Since the code I wrote were **Sage** re-implementations, I have kept those in this repository.

Finally, for all the code files, I left a comment at top indicating whether it references any of the two sources above. We will denote 1. as **PSV** and 2. as **FJSVV**. If there are no comments at top indicaitng this, then that file didn't use any of the two sources.
