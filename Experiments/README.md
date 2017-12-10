# psiTurk implementations of the experiments

the experiments are implemented in HTML & Javascript using the psiTurk framework

they depend on node.js and packages browserify, es6ify

node.js workflow on mac os:
- install node 
```
brew install node
```
(on ubuntu: apt-get install nodejs-legacy)

- install required packages locally:
```
npm install [packagename]
```

- create a bundle.js file, necessary to runthe experiment:
```
node build.js
``` 


psiTurk workflow to test the experiments (linux or mac os only)

- run psiturk in the experiment root directory (psiturk can be installed via pip)
```
psiturk
```

- in the psiturk console, run:
```
server on
debug
```


data
- to pull the data from participants.db (or remote db) into a csv file, first adapt get_data.py as needed (which columns? which table?), then run:
```
python get_data.py
```

(credits to Fabian Dablander)
