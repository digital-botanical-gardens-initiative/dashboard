# dashboard
### *Still in progress*

Set of scripts to render a dashboard that allows to navigate through the lotus dataset.


Two dashboards were created:

## Dash
This dashboard uses the *dash* python library.

To run: 
```bash
cd dash/

conda env create -f environment.yml

python dashboard.py
```

This code isn't being updated anymore.

## Node.js
This dashboard uses *Node.js*.

To run:
```bash
cd node/

npm install package.json

export NODE_TLS_REJECT_UNAUTHORIZED=0

nodemon app
```
