from xml.etree import ElementTree as ET # check that this comes with default python3

""" Currently store everything as a dictionary
"""
def read_xcp(filepath):
    xcp = {
        "SCENE": {},
        "SPECIES": {},
        "ATOMS": {},
        "BONDS": {}
    }
    with open(filepath, 'r', encoding='utf-8') as f:
        data = f.read()
        # load the data here
        root = ET.fromstring(data)
        scene = xcp["SCENE"]
        species = xcp["SPECIES"]
        atoms = xcp["ATOMS"]
        bonds = xcp["BONDS"]
        # for SCENE object throw all objects with a single occurance
        # xml_camera = root.findall("./SCENE/CAMERA")
        # if xml_camera:
        #     scene["CAMERA"] = {}
        #     camera_pos = xml_camera[0].findall("POS")[0].text.split()
        #     camera_rot = xml_camera[0].findall("ROT")[0].text.split()
        #     scene["CAMERA"]["pos"] = tuple([float(x) for x in camera_pos])
        #     scene["CAMERA"]["rot"] = tuple([float(x) for x in camera_rot])
        for spec in root.iter("SPEC"):
            radius = float(spec[0].text)
            sp_id = int(spec.attrib["id"])
            species[sp_id] = {}
            species[sp_id]["label"] = spec.attrib["label"]
            species[sp_id]["scale"] = (radius, radius, radius)
            color = [float(x) for x in spec[1].text.split()]
            species[sp_id]["color"] = tuple(color)
        for atom in root.iter("ATOM"):
            at_id = int(atom.attrib["id"])
            atoms[at_id] = {}
            sp_id = int(atom.attrib["species"])
            atoms[at_id]["spinfo"] = species[sp_id]
            pos = [float(x) for x in atom[0].text.split()]
            atoms[at_id]["pos"] = tuple(pos)
        for bond in root.iter("BOND"):
            bd_id = int(bond.attrib["id"])
            bonds[bd_id] = {}
            bonds[bd_id]["type"] = int(bond.attrib["type"])
            # set keys A, B to reference to the relevant atom entries
            bonds[bd_id]["A"] = atoms[int(bond.attrib["A"])]
            bonds[bd_id]["B"] = atoms[int(bond.attrib["B"])]

    return xcp
    
# QUICK DEBUGGING
if __name__ == "__main__":
    import pprint
    xcp = read_xcp("test.xml")
    pprint.pprint(xcp)