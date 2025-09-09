import xml.etree.ElementTree as ET
import json

input_xml = "/u/home/t/tchhabri/project-kappel/clinvar/ClinVarVCVRelease_00-latest.xml"
output_json = "/u/home/t/tchhabri/project-kappel/clinvar/full_output_latest_august_2025.json"

def ET_to_dict(elem):
    d = {}
    for k, v in elem.attrib.items():
        d[f"@{k}"] = v
    text = (elem.text or "").strip()
    if text and not list(elem):
        d["#text"] = text
    children = {}
    for child in elem:
        child_dict = ET_to_dict(child)
        if child.tag not in children:
            children[child.tag] = []
        children[child.tag].append(child_dict)
    for tag, vals in children.items():
        if len(vals) == 1:
            d[tag] = vals[0]
        else:
            d[tag] = vals
    return d

def flatten_dict(d, parent_key='', sep='_'):
    items = []
    if isinstance(d, list):
        for i, item in enumerate(d):
            new_key = f"{parent_key}{sep}{i}" if parent_key else str(i)
            items.extend(flatten_dict(item, new_key, sep=sep).items())
    elif isinstance(d, dict):
        for k, v in d.items():
            new_key = f"{parent_key}{sep}{k}" if parent_key else k
            if isinstance(v, (dict, list)):
                items.extend(flatten_dict(v, new_key, sep=sep).items())
            else:
                items.append((new_key, v))
    else:
        items.append((parent_key, d))
    return dict(items)

def stream_parse_clinvar(xml_path, json_path):
    with open(json_path, "w", encoding="utf-8") as out_f:
        out_f.write("[\n")
        first = True
        context = ET.iterparse(xml_path, events=("end",))
        for event, elem in context:
            if elem.tag == "VariationArchive":
                nested_dict = ET_to_dict(elem)
                flat_dict = flatten_dict(nested_dict)
                if not first:
                    out_f.write(",\n")
                else:
                    first = False
                json.dump(flat_dict, out_f, indent=2)
                elem.clear()
        out_f.write("\n]")

stream_parse_clinvar(input_xml, output_json)

