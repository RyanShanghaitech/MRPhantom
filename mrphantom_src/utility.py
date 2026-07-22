import xml.etree.ElementTree as ET

def ele2dict(element):
    result = {
        key: int(value) if key=="enum" else float(value)
        for key, value in element.attrib.items()
    }
    for child in element:
        result[child.tag] = ele2dict(child)
    return result

xml2dict = lambda fXml: ele2dict(ET.parse(fXml).getroot())