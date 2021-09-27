
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

var realConsoleLog = console.log;
console.log = function () {
    var message = [].join.call(arguments, " ");
    postMessage(message);
    realConsoleLog.apply(console, arguments);
};

onmessage = function(e) {
    
    var filenames = []
    for(var i =0; i<e.data[0].length; i++){
        filenames.push(e.data[0][i].name);
    }
    

    // print some log info
    console.log("filenames: ", filenames);
    console.log("e.data[1]: ", e.data[1]);
    console.log("e.data[2]: ", e.data[2]);
    console.log("e.data[3]: ", e.data[3]);
    console.log("e.data[4]: ", e.data[4]);
    console.log("e.data[5]: ", e.data[5]);
    console.log("e.data[6]: ", e.data[6]);
    console.log("e.data[7]: ", e.data[7]);
    console.log("e.data[8]: ", e.data[8]);
    console.log("e.data[9]: ", e.data[9]);
    console.log("e.data[10]: ", e.data[10]);

    // mount the input files
    console.log("Loading file...");
    if (!FS.analyzePath('/work').exists){
        FS.mkdir('/work');
    }
    FS.mount(WORKERFS, { files: e.data[0] }, '/work');


    if (e.data[10]==='snp'){
        // run pairsnp
        var retVector = Module.pairsnp('/work/' + filenames[0], 
        1, parseInt(e.data[8]), parseInt(e.data[9]));
    } else if (e.data[10]==='accessory') {
        // run pairwise binary distance
        var retVector = Module.pairgene('/work/' + filenames[0], 
            1, parseInt(e.data[8]), parseInt(e.data[9]));
    } else if (e.data[10]==='sketch') {
        // run pairwise sketch distance
        console.log('HEREHERHEHRHEHR!!')
        var vec = new Module.StringList();
        for(var i =0; i<e.data[0].length; i++){
            vec.push_back('/work/'  + filenames[i]);
            // filenames[i] = '/work/'  + filenames[i];
            // console.log("filenames[i]: ", filenames[i]);
        }
        var retVector = Module.ppsketch(vec, 
            1, -1, -1, 0);
    }

    // run SCE
    var wtsneVector = [...vector_values(Module.wtsne(retVector.rows,
        retVector.cols, 
        retVector.distances,
        retVector.seq_names.size(),
        parseFloat(e.data[1]),
        parseInt(e.data[2]),
        parseInt(e.data[3]),
        parseInt(e.data[4]),
        parseFloat(e.data[5]),
        parseInt(e.data[6]),
        parseInt(e.data[7])))];

    postMessage(wtsneVector);

    FS.unmount('/work');
    
}

self.importScripts('pathosce.js');