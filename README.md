# CDCgov GitHub Organization Open Source Project 



**General disclaimer** This repository was created for use by CDC programs to collaborate on public health related projects in support of the [CDC mission](https://www.cdc.gov/about/organization/mission.htm).  GitHub is not hosted by the CDC, but is a third party website used by CDC and its partners to share information and collaborate on software. CDC use of GitHub does not imply an endorsement of any one particular service, product, or enterprise. 



## Overview

This bioinformatics pipeline identifies and summarizes amplicon sequences produced by the multiplex PCR amplicon sequencing (MPAS) described by Hojgaard et. al (2020). The MPAS assay was originally designed to detect microbial genera that contain known human pathogens found in _Ixodes_ ticks (_Borrelia_ spp. _Babesia_ spp., _Anaplasma_ spp. _Ehrlichia_ spp.). This pipeline is implemented in Nextflow and performs QC analysis, primer trimming, error correction, alignment to reference sequences, NCBI BLAST of unaligned sequences, phylogenetic tree creation of aligned sequences, and report generation.

To run this pipeline, see [pipeline instructions](Pipeline_instructions.md)


This pipeline was developed by Mark Stenglein from Colorado State University in collaboration with Lynn Osikowicz, Sarah Maes, Andrias Hojgaard, and Becky Eisen from the CDC's Division of Vector-Borne Diseases. 

  
## Public Domain Standard Notice
This repository constitutes a work of the United States Government and is not
subject to domestic copyright protection under 17 USC § 105. This repository is in
the public domain within the United States, and copyright and related rights in
the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
All contributions to this repository will be released under the CC0 dedication. By
submitting a pull request you are agreeing to comply with this waiver of
copyright interest.

## License Standard Notice
The repository utilizes code licensed under the terms of the Apache Software
License and therefore is licensed under ASL v2 or later.

This source code in this repository is free: you can redistribute it and/or modify it under
the terms of the Apache Software License version 2, or (at your option) any
later version.

This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the Apache Software License for more details.

You should have received a copy of the Apache Software License along with this
program. If not, see http://www.apache.org/licenses/LICENSE-2.0.html

The source code forked from other open source projects will inherit its license.

## Privacy Standard Notice
This repository contains only non-sensitive, publicly available data and
information. All material and community participation is covered by the
[Disclaimer](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md)
and [Code of Conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).
For more information about CDC's privacy policy, please visit [http://www.cdc.gov/other/privacy.html](https://www.cdc.gov/other/privacy.html).

## Contributing Standard Notice
Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo)
and submitting a pull request. (If you are new to GitHub, you might start with a
[basic tutorial](https://help.github.com/articles/set-up-git).) By contributing
to this project, you grant a world-wide, royalty-free, perpetual, irrevocable,
non-exclusive, transferable license to all users under the terms of the
[Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or
later.

All comments, messages, pull requests, and other submissions received through
CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

## Records Management Standard Notice
This repository is not a source of government records, but is a copy to increase
collaboration and collaborative potential. All government records will be
published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices
Please refer to [CDC's Template Repository](https://github.com/CDCgov/template)
for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/master/CONTRIBUTING.md),
[public domain notices and disclaimers](https://github.com/CDCgov/template/blob/master/DISCLAIMER.md),
and [code of conduct](https://github.com/CDCgov/template/blob/master/code-of-conduct.md).

## Pipeline References  
Hojgaard, A., Osikowicz, L. M., Eisen, L., & Eisen, R. J. (2020). Evaluation of a novel multiplex PCR amplicon sequencing assay for detection of human pathogens in Ixodes ticks. Ticks and tick-borne diseases, 11(6), 101504. doi:10.1016/j.ttbdis.2020.101504.  [PMID: 32993925](https://pubmed.ncbi.nlm.nih.gov/32993925/).

Osikowicz LM, Hojgaard A, Maes S, Eisen RJ, Stenglein MD. A bioinformatics pipeline for a tick pathogen surveillance multiplex amplicon sequencing assay.  Ticks Tick Borne Dis. 2023 Sep;14(5):102207. doi: 10.1016/j.ttbdis.2023.102207.  Epub 2023 May 27. [PMID: 37247570](https://pubmed.ncbi.nlm.nih.gov/37247570/).
