initRDKit = async function initRDKit(smile, id) {
    await Module.initRDKit();
    const mol = Module.get_mol(smile);
    const svg = mol.get_svg();
    document.getElementById(id).innerHTML = svg;
    mol.delete();
  };

module.exports = {
    initRDKit: initRDKit
};