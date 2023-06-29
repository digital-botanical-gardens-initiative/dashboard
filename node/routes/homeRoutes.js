// homeRoutes.js
const express = require('express');
const router = express.Router();
const db = require('./db');

router.use(express.urlencoded({ extended: true }));

router.get('/', async (req, res) => {

    const plants = await db.query(`
        SELECT DISTINCT organism_taxonomy_08genus
        FROM data
        WHERE organism_taxonomy_02kingdom = 'Archaeplastida'`);


    
/*     const resultMol = await db.query(`
        SELECT 
            structure_taxonomy_npclassifier_01pathway as pathway,
            structure_taxonomy_npclassifier_02superclass as superclass,
            structure_taxonomy_npclassifier_03class as class
        FROM data
        GROUP BY pathway, superclass, class
    `);

    const dataMol = resultMol.rows;

    const treeMol = { name: "Root", children: [] };
    dataMol.forEach(row => {
        let currentLevel = treeMol.children;
        [row.pathway, row.superclass, row.class].forEach(level => {
            let existingPath = currentLevel.find(d => d.name === level);
            if (existingPath) {
                currentLevel = existingPath.children;
            } else {
                const newPath = { name: level, children: [] };
                currentLevel.push(newPath);
                currentLevel = newPath.children;
            }
        });
    });


    const resultOrg = await db.query(`
        SELECT 
        organism_taxonomy_01domain as domain,
        organism_taxonomy_02kingdom as kingdom,
        organism_taxonomy_03phylum as phylum,
        organism_taxonomy_04class as class,
        organism_taxonomy_06family as family,
        organism_taxonomy_07tribe as tribe,
        organism_taxonomy_08genus as genus
        FROM data
        GROUP BY domain, kingdom, phylum, class, family, tribe, genus
    `);

    const dataOrg = resultOrg.rows;

    const treeOrg = { name: "Root", children: [] };
    dataOrg.forEach(row => {
        let currentLevel = treeOrg.children;
        [row.domain, row.kingdom, row.phylum, row.class, row.family,
        row.tribe, row.genus].forEach(level => {
            let existingPath = currentLevel.find(d => d.name === level);
            if (existingPath) {
                currentLevel = existingPath.children;
            } else {
                const newPath = { name: level, children: [] };
                currentLevel.push(newPath);
                currentLevel = newPath.children;
            }
        });
    });

    res.render('home', { treeMol, treeOrg }); */

    var result = plants.rowCount;

    res.render('home', {result});
});


module.exports = router;

