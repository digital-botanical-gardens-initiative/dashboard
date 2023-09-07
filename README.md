# DBGI dashboard
### *Still in progress*

Set of scripts to render a dashboard that allows to navigate through the LOTUS and DBGI dataset.


## Abstract of the project
The Digital Botanical Gardens Initiative (DBGI) embarks on an innovative journey to curate, manage, and disseminate digital data from living botanical collections, with an emphasis on mass spectrometric evaluations of chemodiversity. Using semantic web technology, this data is linked with relevant metadata, propelling ecosystem research and guiding biodiversity conservation efforts. Central to the success of DBGI is the creation of an interactive platform for both humans and machines to assimilate this knowledge. This report outlines our efforts to design the prototype of a data visualization portal intended to evolve into the DBGI dashboard. Starting with a Plotly Dash application, the project transitioned to a Node.js application leveraging Javascript, HTML, and CSS for enhanced customization. This provides a basis for future improvements, some of which are proposed in the report.

## 
**Two dashboards were created for the project:**

## Dash
This dashboard uses the *dash* python library and only navigates through the LOTUS dataset.

To run: 
```bash
cd dash/

conda env create -f environment.yml

python dashboard.py
```

**This code isn't being updated anymore.**

## Node.js
This dashboard uses *Node.js* and navigates through the LOTUS and the DBGI dataset.

### Installation:

1. Clone repository

2. Install packages

```bash
cd node/

npm install package.json
```

3. Run the app
```bash
nodemon app
```

### Features

The dashboard is organized as a multi-page application, enabling users to explore both databases via various search methods.

The top of each page features a navigation bar, ensuring accessibility to each section from any page within the app. The navigation bar contains the title, which when clicked, redirects the user to the home page. Two buttons are present, one also directing the user to the home page and the other to the download page.

Furthermore, there are two dropdown menus. The first, titled *Explore*, allows users to navigate the databases in different ways, such as *Explore Whole Dataset*, *Explore Molecular Diversity*, and *Explore Organisms Diversity*.

The second dropdown menu, *Information*, provides access to essential documentation and background information about the DBGI. It also houses a link to the DBGI organization's GitHub page, enabling direct access to the project's codebase and updates.

#### Home
The Home page serves as an introductory interface for the user, from where they can navigate to various sections of the application.

Additionally, it provides insightful project statistics, such as the proportion of various taxonomic categories - including phylum, class, and species - that have already been profiled within the project. This feature offers users a snapshot view of the current state of the project's data collection and profiling progress.

The values for the total number of taxon known were taken from the *Catalogue of Life* website.

#### Explore whole dataset
This section gives the user the ability to access and navigate through the LOTUS and DBGI databases using a variety of methods. The primary methods include a text-based search where the user can hunt for specific words within the databases, a structural search that allows the user to sketch a molecule's structure and look for it, and direct SPARQL queries. It should be noted, however, that direct SPARQL queries can only be used with the DBGI data because the LOTUS database is housed in a relational database that relies on SQL.

The user can choose which dataset to browse by selecting the datasource in both text- and structure-based search. When browsing the DBGI data, the back-end queries are SPARQL queries and when browsing the LOTUS data, the back-end queries are SQL queries.

The results from both search methods can be visualized in two distinct formats based on user preference: tabular form or treemap representation. The treemaps are generated using the Plotly library for JavaScript. The organization of the treemaps mirrors the phylogeny of the organisms where the molecule or substructure can be identified. A colorblind-accessible color scheme is employed for this visualization. Owing to Plotly's interactive features, users can engage with the graph, zooming in on sections that pique their interest.

The output of the search operation is constrained to a user-defined quantity of rows.

The next sections go through the specificity of each search method.

**Text-based**      

The text-based search operation conducts a case-insensitive scan for the text input within a user-selected column. The resulting output is contingent on the specified row limit and the chosen display format.

**Structure-based**

In the structure-based search page, the user can draw a molecule and then search for hits using different methods. 
The JSME object provides a variety of shortcuts to facilitate user interactions when drawing molecules. Detailed instructions for these can be accessed on the \href{https://jsme-editor.github.io/help.html}{JSME help page}.

One of these shortcuts allows users to insert the SMILES notation of the molecule they want to represent, eliminating the need for manual drawing.

##### Exact Match

The search for an exact match can be accomplished using either the InChI or SMILES representation of the drawn molecule, each offering distinct advantages:

For optimal searching accuracy, InChI is preferable as it minimizes the chance of ambiguity during searches.

Based on the selected representation, the dashboard will query the designated database to retrieve all relevant matches where the InChI or the SMILES align with the depicted structure.

##### Substructure Search

The structure-based search feature empowers users to query for molecules that contain a specific drawn substructure.

In practice, the user drafts a unique structural design, which the application then uses as a query to search the database for molecules incorporating this substructure.

##### Similarity Search

The Jaccard-Tanimoto coefficient is utilized in conducting the similarity search. Users can choose the percentage of similarity they want between the drawn structure and the compounds queried.

**SPARQL**                      

The concluding segment of the *Explore Whole Dataset* division involves the SPARQL search.

This section empowers users to directly interface with the DBGI Knowledge Graph by crafting and executing their own SPARQL queries.

#### Explore Molecules
This segment of the application allows users to query the LOTUS dataset for specific molecules via their names. Additionally, it showcases the 12 most frequently encountered molecules within the LOTUS database, thus highlighting key molecules of interest.

The molecules name displayed are links to a page described in next section, that display informations about the molecule.

#### Molecule Page
Upon choosing a molecule from either the *Explore Molecules* page or the *Explore* table results, the user is redirected to a page offering detailed information about the molecule.

Adjacent to the name of the molecule, the Wikidata ID is displayed, which serves as a direct link to the Wikidata entity of the respective molecule. Following this, the 2D and 3D structures of the molecule are presented. Next, the chemical taxonomy (obtained using Classyfire [1]) of the molecule is listed.

Lastly, the page provides a list of organisms in which this molecule has been identified in the LOTUS database. Furthermore, the names of these organisms serve as links to respective organism pages, offering an in-depth view into each organism.

#### Explore Organism
Similar to the *Explore Molecules* section, this segment provides the user the ability to search for a specific organism within the LOTUS database.

In order to aid the user and provide a starting point, a list of the top 12 most prevalent organisms is also displayed.

Importantly, each result acts as a direct link to a dedicated page, described in the next section, for the specific organism, offering a deeper exploration of each individual organism.

#### Organism page
This page presents information pertaining to organisms catalogued in the LOTUS database, with each page dedicated to a specific organism. Access to these pages can be achieved through either the *Explore Organism* page or the *Explore* table results.

Alongside the organism's name, its Wikidata ID is displayed, which serves as a direct link to the corresponding Wikidata entity page.

Furthermore, the Wikidata image of the organism is displayed below, accompanied by its phylogenetic information.

Lastly, a list of all molecules identified within this organism, as recorded in the LOTUS database, is displayed. Notably, each molecule name serves as a hyperlink to a dedicated page for that specific molecule, thereby facilitating further exploration.

#### Download
This segment provides users with the capability to download the LOTUS dataset in either CSV or JSON format. Regrettably, this feature is currently non-operational due to an issue that causes the application to crash upon attempted download. Work is underway to rectify this problem.


## References
[1] Djoumbou Feunang, Y., Eisner, R., Knox, C. et al. ClassyFire: automated chemical classification with a comprehensive, computable taxonomy. J Cheminform 8, 61 (2016). https://doi.org/10.1186/s13321-016-0174-y
