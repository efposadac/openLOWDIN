from xml.dom import minidom

xmldoc = minidom.parse('ee.xml')

itemlist = xmldoc.getElementsByTagName('Interaction')

for item in itemlist:
    print("O-" + item.attributes['name'].value)
    print("# Add commentary")
    print(item.attributes['size'].value)
    gclist = item.getElementsByTagName('GaussianComponent')
    for i, gc in enumerate(gclist):
        print(i, int(float(gc.attributes['angularMoment'].value)))
        print(gc.attributes['exponent'].value, gc.attributes['factor'].value)
        print(gc.attributes['x'].value, gc.attributes[
              'y'].value, gc.attributes['z'].value)

    print("")
