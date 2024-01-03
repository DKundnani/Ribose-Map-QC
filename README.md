
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
    <li><a href="#Usage">Usage</a></li>
      <ul>
        <li><a href="###Getting QC report">Getting QC report</a></li>
       <li><a href="###Getting-counts-using-bed-file">Getting-counts-using-bed-file</a></li>
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
conda env create --file /Ribose-Map-QC/env.yml
```

## Getting QC report
### Configure run
```bash
vim MMremoval_configure_run.sh
#Change the variables for you run as per mentions in the bash configure file
```

### Filtration
```bash
bash MMremoval_configure_run.sh

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

Deepali L. Kundnani- [![LinkedIn][linkedin-shield]][linkedin-url] - [deepali.kundnani@gmail.com](mailto::deepali.kundnani@gmail.com)

Project Link: [https://github.com/your_username/repo_name](https://github.com/your_username/repo_name)

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

Use this space to list resources you find helpful and would like to give credit to. I've included a few of my favorites to kick things off!

* [Choose an Open Source License](https://choosealicense.com)
* [GitHub Emoji Cheat Sheet](https://www.webpagefx.com/tools/emoji-cheat-sheet)
* [Malven's Flexbox Cheatsheet](https://flexbox.malven.co/)
* [Malven's Grid Cheatsheet](https://grid.malven.co/)
* [Img Shields](https://shields.io)
* [GitHub Pages](https://pages.github.com)
* [Font Awesome](https://fontawesome.com)
* [React Icons](https://react-icons.github.io/react-icons/search)

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