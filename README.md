### JS Export of Simplest SIR Model with Hyperparameter Tuning: https://silico.app/@leonfan/sir-for-hyperparameter-tuning?s=dp6qsWdmQbyjpyhoBwieJA 

## How to run the model:
- Navigate into COVID-19 directory
- node model2.js (make sure node is installed on your computer)

## Helpful functions within model2.js
- Optimize (takes in 2 lists, uses time series to tune second parameter)
- Optimize_simple (takes in 2 lists, uses flat values to tune)
- loss (calculates loss via least squares)
- test (run the model to access simulation data)
- export_to_json (results in output.json file that contains hyperparameter values)

## See all annotations within model2.js

## Things to note:
- Within terminal, you should see console.log("complete") for successful runs

