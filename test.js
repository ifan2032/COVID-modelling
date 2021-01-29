var fs = require("fs");

data = {arr: [1,2]};
fs.writeFile ("output.json", JSON.stringify(data), (err) => {
    if (err) throw err;
    console.log("complete");
})
