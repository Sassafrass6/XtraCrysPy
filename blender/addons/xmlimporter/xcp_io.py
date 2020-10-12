from xml.etree import ElementTree as ET # check that this comes with default python3

""" Currently store everything as a dictionary
"""

def get_tag(root, tag):
    try:
        return root.findall(tag)[0]
    except ValueError:
        return None
    except IndexError:
        return None

def get_attrib(root, tag, att, typ):
    try:
        if typ == str:
            return root.findall(tag)[0].attrib[att]
        elif typ == int:
            return int(root.findall(tag)[0].attrib[att])
        elif typ == float:
            return float(root.findall(tag)[0].attrib[att])
        elif typ == tuple:
            return tuple([float(x) for x in root.findall(tag)[0].attrib[att].split()])
    except ValueError:
        return None
    except IndexError:
        return None
    except KeyError:
        return None
    return None
    
def get_text(root, tag, typ):
    try:
        if typ == str:
            return root.findall(tag)[0].text
        elif typ == int:
            return int(root.findall(tag)[0].text)
        elif typ == float:
            return float(root.findall(tag)[0].text)
        elif typ == tuple:
            return tuple([float(x) for x in root.findall(tag)[0].text.split()])
    except ValueError:
        return None
    except IndexError:
        return None
    return None

def read_xcp(filepath):
    xcp = {
        "SCENE": {},
        "SPECIES": {},
        "ATOMS": {},
        "BONDS": {},
        "FRAME": {}
    }
    with open(filepath, 'r', encoding='utf-8') as f:
        data = f.read()
        # load the data here
        root = ET.fromstring(data)
        scene = xcp["SCENE"]
        species = xcp["SPECIES"]
        atoms = xcp["ATOMS"]
        bonds = xcp["BONDS"]
        frame = xcp["FRAME"]
        # for SCENE object throw all objects with a single occurance
        if get_tag(root, './SCENE/CAMERA'):
            scene["CAMERA"] = {}
            scene["CAMERA"]["position"] = get_text(root, './SCENE/CAMERA/POSITION', tuple)
            scene["CAMERA"]["rotation"] = get_text(root, './SCENE/CAMERA/ROTATION', tuple)
        for spec in root.iter("SPECIES"):
            sp_id = get_attrib(spec, '.', 'id', int)
            species[sp_id] = {}
            species[sp_id]["label"] = get_attrib(spec, '.', 'label', str)
            species[sp_id]["scale"] = get_text(spec, './RADIUS', float)
            species[sp_id]["color"] = get_text(spec, './COLOR', tuple)
        for atom in root.iter("ATOM"):
            at_id = get_attrib(atom, '.', 'id', int)
            atoms[at_id] = {}
            sp_id = get_attrib(atom, '.', 'species', int)
            atoms[at_id]["spinfo"] = species[sp_id]
            atoms[at_id]["position"] = get_text(atom, '.', tuple)
        for bond in root.iter("BOND"):
            bd_id = get_attrib(bond, '.', 'id', int)
            bonds[bd_id] = {}
            bonds[bd_id]["type"] = get_attrib(bond, '.', 'type', int)
            # set keys A, B to reference to the relevant atom entries
            bonds[bd_id]["A"] = atoms[get_attrib(bond, '.', 'A', int)]
            bonds[bd_id]["B"] = atoms[get_attrib(bond, '.', 'B', int)]
        xframe = get_tag(root, './FRAME')
        if xframe:
            # expect vertices and edges
            frame["VERTICES"] = {}
            frame["EDGES"] = {}
            for vertex in xframe.iter("VERTEX"):
                vx_id = get_attrib(vertex, '.', 'id', int)
                frame["VERTICES"][vx_id] = { "position": get_text(vertex, '.', tuple) }
            for edge in xframe.iter("EDGE"):
                eg_id = get_attrib(edge, '.', 'id', int)
                frame["EDGES"][eg_id] = { 
                    "A": frame["VERTICES"][get_attrib(edge, '.', 'A', int)],
                    "B": frame["VERTICES"][get_attrib(edge, '.', 'B', int)]
                }

    return xcp
    
# QUICK DEBUGGING
if __name__ == "__main__":
    import pprint
    xcp = read_xcp("test.xml")
    pprint.pprint(xcp)