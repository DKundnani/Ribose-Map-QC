
<h1 align="center">Ribose-Map-QC</h1>
To get quality control of files generated by Ribose-Map. 
Also supported by percentage and composition of rNMPs from ribose-seq output 

<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--
-->
[![Commits][Commits-shield]][Commits-url]
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Website][website-shield]][website-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]

<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li><a href="##Installation">Installation</a></li>
      <ul>
        <li><a href="###getting-the-code">Getting the code</a></li>
        <li><a href="###Creating-the-enviroment-with-required-dependencies">Creating the enviroment with required dependencies</a></li>
      </ul>
    </li>
    <li><a href="##Usage">Usage</a></li>
      <ul>
        <li><a href="###Getting QC report">Getting QC report</a></li>
        <li><a href="###Getting-counts-using-bed-file">Getting counts using bed file</a></li>
       <li><a href="###Visualization">Visualization</a></li>
      </ul>
  </ol>
</details>


## Installation

### Getting the code
The development version from [GitHub](https://github.com/) with:

```sh
git clone https://github.com/DKundnani/Ribose-Map-QC.git
```

### Creating the enviroment with required dependencies

```sh
conda env create --name RibosemapQC_env --file /Ribose-Map-QC/env.yml
```

### Additional Dependencies

* Ribosemap generated output
* Trimmed reads used in Ribosemap
* Reference genome of the organism being used (also used in Ribosemap)

## Usage
### Getting QC report
```bash
cd $outputdir #This would be the results directory generated by Ribose-Map. Please make sure you have generate alignment, coordinate and composition modules to get the results and the names to trimmed files match names of folders as output of Ribosemap. files will contain the library name as first column and sample name as second column(example present in repository.
python3 path/to/Ribose-Map-QC/quality.py -f path/to/files -t path/to/trimmed_reads -r path/to/reference/sacCer2/sacCer2.fa #Generates QC.tsv as output 
```

### Getting counts using bed file
```bash
Place all the bed files in a folder e.g. path/to/bed
generate_QC_frombed.py f path/to/files -r path/to/reference/sacCer2/sacCer2.fa #Generates QC.tsv as output
```

### Visualization
```bash
conda activate RibosemapQC_env #activating Enviroment
Rscript path/to/Ribose-Map-QC/qualviz.R -q QC.tsv -l files #Using QC.tsv generated from either quality.py or qual_short.py
Rscript path/to/Ribose-Map-QC/counts_viz.R -c QC.tsv -l files #Using QC.tsv generated from generate_QC_frombed.py
```





<!-- CONTRIBUTING -->
## Contributing

Contributions are what make the open source community such an amazing place to learn, inspire, and create. Any contributions you make are **greatly appreciated**.

If you have a suggestion that would make this better, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".
Don't forget to give the project a star! Thanks again!

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/AmazingFeature`)
3. Commit your Changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the Branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

Deepali L. Kundnani - [deepali.kundnani@gmail.com](mailto::deepali.kundnani@gmail.com)    [![LinkedIn][linkedin-shield]][linkedin-url] 

#Project Link: [https://github.com/your_username/repo_name](https://github.com/your_username/repo_name)

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- ACKNOWLEDGMENTS -->
## Citations

Use this space to list resources you find helpful and would like to give credit to. I've included a few of my favorites to kick things off!

* <b>Distinct features of ribonucleotides within genomic DNA in Aicardi-Goutières syndrome (AGS)-ortholog mutants of <i>Saccharomyces cerevisiae</i> </b>
Deepali L. Kundnani, Taehwan Yang, Alli L. Gombolay, Kuntal Mukherjee, Gary Newnam, Chance Meers, Zeel H. Mehta, Celine Mouawad, Francesca Storici
bioRxiv 2023.10.02.560505; doi: https://doi.org/10.1101/2023.10.02.560505


<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/DKundnani/Ribose-Map-QC?style=for-the-badge
[contributors-url]: https://github.com/DKundnani/Ribose-Map-QC/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/DKundnani/Ribose-Map-QC?style=for-the-badge
[forks-url]: https://github.com/DKundnani/Ribose-Map-QC/forks
[stars-shield]: https://img.shields.io/github/stars/DKundnani/Ribose-Map-QC?style=for-the-badge
[stars-url]: https://github.com/DKundnani/Ribose-Map-QC/stargazers
[issues-shield]: https://img.shields.io/github/issues/DKundnani/Ribose-Map-QC?style=for-the-badge
[issues-url]: https://github.com/DKundnani/Ribose-Map-QC/issues
[license-shield]: https://img.shields.io/github/license/DKundnani/Ribose-Map-QC?style=for-the-badge
[license-url]: https://github.com/DKundnani/Ribose-Map-QC/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/deepalik
[product-screenshot]: images/screenshot.png
[commits-url]: https://github.com/DKundnani/Ribose-Map-QC/pulse
[commits-shield]: https://img.shields.io/github/commit-activity/t/DKundnani/Ribose-Map-QC?style=for-the-badge
[website-shield]: https://img.shields.io/website?url=http%3A%2F%2Fdkundnani.bio%2F&style=for-the-badge
[website-url]:http://dkundnani.bio/ 
