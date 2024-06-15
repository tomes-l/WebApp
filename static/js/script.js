// Function to handle SMILES form submission
const formSMILES = document.getElementById('form-smiles')

formSMILES.addEventListener('submit', event => {
  event.preventDefault(); // Prevent default submission

  const formData = new FormData(formSMILES);  // Create dict of keys/values
  const SMILESinput = Object.fromEntries(formData); //Formate the dict in Js object

  fetch('/submit-form-smiles', {
    method: 'POST',
    headers: {
        'Content-Type' : 'application/json'
    },
    body: JSON.stringify(SMILESinput)
  })
  .then(res => {
    return res.json()
  })
  .then(data => {
    console.log(data)
    if (data.img) {
      const imgRes = data.img;
      const img = document.createElement('img');
      img.src = 'data:image/png;base64,' + imgRes;
      document.getElementById('representation-from-SMILES').appendChild(img);
    }

    if (data.mol) {
      const MolRes = document.createElement('p');
      MolRes.innerHTML = `Mol: ${data.mol}`;
      document.getElementById('results-smiles').appendChild(MolRes);
    }

    if (data.mformula) {
      const MFRes = document.createElement('p');
      MFRes.innerHTML = `Molecular formula: ${data.mformula}`;
      document.getElementById('results-smiles').appendChild(MFRes);
    }

    if (data.mweight) {
      const MWres = document.createElement('p');
      MWres.innerHTML = `Molecular weight: ${data.mweight}`;
      document.getElementById('results-smiles').appendChild(MWres);
    }
  })
  .catch(error => {
  const ErrorSMILES = document.createElement('p');
  ErrorSMILES.innerHTML = `There was an error processing your request`;
  document.getElementById('results-mol').appendChild(ErrorSMILES);
  })
});

//////////////////////////////////////
////////////Molfile Part//////////////
//////////////////////////////////////

// Function to handle MOL file
function readFileAsString() {
  return new Promise((resolve, reject) => {
      const [file] = document.getElementById('Molfile-input').files;
      const reader = new FileReader();

      reader.addEventListener("load", () => {
          resolve(reader.result); // Solve the promise with the contents of the file
      });

      reader.addEventListener("error", () => {
          reject(reader.error); // Reject promise in case of errors
      });

      if (file) {
          reader.readAsText(file); // Read the file if select
      } else {
          reject("No file select"); // Reject promise if no file select
      }
  });
}


// Event listener for Mol form when submit click
const formMol = document.getElementById('form-mol')

formMol.addEventListener('submit', async (event) => {
  event.preventDefault();
  
  const file_str = await readFileAsString(); 
  console.log(file_str)
  const formData = new FormData(formMol);
  const Molinput = Object.fromEntries(formData); //Formate the dict in Js object
  console.log(Molinput)

  let dictMol = {};
  

  // Go through the form elements
  // The goal here is to stock the value null to text input if they are empty
  // For the back-end python to well operate on the inputs data
  for (let [key, value] of Object.entries(Molinput)) {
    const inputElement = formMol.querySelector(`[name="${key}"]`);
    if (inputElement && inputElement.type === 'text' && value === '') {
        dictMol[key] = null ; // If the input is empty, store null
    } else {
        dictMol[key] = value; // Otherwise, stores the value of the input
    }
  }

  dictMol['file'] = file_str;
  console.log(JSON.stringify(dictMol))

  fetch('/submit-form-mol', {
    method: 'POST',
    headers: {
        'Content-Type' : 'application/json'
    },
    body: JSON.stringify(dictMol)
  })
  .then(res => {
    return res.json()
  })
  .then(data => {
    console.log(data)
    if (data.img) {
      // Print the molecular representation
      const imgRes = data.img;
      const img = document.createElement('img');
      img.src = 'data:image/png;base64,' + imgRes;
      document.getElementById('representation-from-mol').appendChild(img);
    }

    if (data.smiles) {
      const smilesRes = document.createElement('p');
      smilesRes.innerHTML = `SMILES: ${data.smiles}`;
      document.getElementById('results-mol').appendChild(smilesRes);
    }

    if (data.count_atom) {
      const CountAtomRes = document.createElement('p');
      CountAtomRes.innerHTML = `Count Atom: ${data.count_atom}`;
      document.getElementById('results-mol').appendChild(CountAtomRes);
    }

    if (data.bond_distance) {
      const BondDistanceRes = document.createElement('p');
      BondDistanceRes.innerHTML = `Topological distance: ${data.bond_distance}`;
      document.getElementById('results-mol').appendChild(BondDistanceRes);
    }

    if (data.room_distance) {
      const RoomDistanceRes = document.createElement('p');
      RoomDistanceRes.innerHTML = `3D distance: ${data.room_distance}`;
      document.getElementById('results-mol').appendChild(RoomDistanceRes);
    }

    if (data.atom_neighbours) {
      const AtomNeighboursRes = document.createElement('p');
      AtomNeighboursRes.innerHTML = `Atom neighbours: ${data.atom_neighbours}`;
      document.getElementById('results-mol').appendChild(AtomNeighboursRes);
    }

    if (data.ring) {
      const RingRes = document.createElement('p');
      RingRes.innerHTML = `Ring: ${data.ring}`;
      document.getElementById('results-mol').appendChild(RingRes);
    }

    if (data.error) {
      const ErrorRes = document.createElement('p');
      ErrorRes.innerHTML = `Error: ${data.error}`;
      document.getElementById('results-mol').appendChild(ErrorRes);
    }
  })
  .catch(error => {
    const ErrorMol = document.createElement('p');
    ErrorMol.innerHTML = `There was an error processing your request`;
    document.getElementById('results-mol').appendChild(ErrorMol);
  })
});

