
function* vector_values(vector) {
    var matrix = []
    for (let k=0; k<2; k++){
        matrix[k] = [];
        for (let i = k; i < vector.size(); i+=2)
            matrix[k].push(vector.get(i));
    }
    yield(matrix);
    vector.delete();
}


onmessage = function(e) {
    const f = e.data[0];

    FS.mkdir('/work');
    FS.mount(WORKERFS, { files: [f] }, '/work');

    var retVector = Module.pairsnp('/work/' + f.name, 1, 10, 11);
    // vector size
    var vectorSize = retVector.rows.size();

    // retrieve value from the vector
    console.log("file: ", f.name);
    console.log("e.data[1]: ", e.data[1]);
    console.log("e.data[2]: ", e.data[2]);
    console.log("e.data[3]: ", e.data[3]);
    console.log("e.data[4]: ", e.data[4]);
    console.log("e.data[5]: ", e.data[5]);
    console.log("e.data[6]: ", e.data[6]);
    console.log("e.data[7]: ", e.data[7]);

    var wtsneVector = [...vector_values(Module.wtsne(retVector.rows,
        retVector.cols, 
        retVector.distances,
        parseFloat(e.data[1]),
        parseInt(e.data[2]),
        parseInt(e.data[3]),
        parseInt(e.data[4]),
        parseFloat(e.data[5]),
        parseInt(e.data[6]),
        parseInt(e.data[7])))];

    postMessage(wtsneVector);
    
}

self.importScripts('/js/pathosce.js');