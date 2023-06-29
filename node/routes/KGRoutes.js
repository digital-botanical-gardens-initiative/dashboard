// downloadRoutes.js
const express = require('express');
const router = express.Router();
const db = require('./db');

router.use(express.urlencoded({ extended: true }));

router.get('/KG', (req, res) => {
    res.render('KG', {title: 'KG'});
});

module.exports = router;

