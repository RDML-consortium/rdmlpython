#!/usr/bin/python

import sys
import os
import re
import datetime
import zipfile
import argparse
from lxml import etree as ET


class RdmlError(Exception):
    """Basic exception for errors raised by the RDML-Python library"""
    def __init__(self, message):
        Exception.__init__(self, message)
    pass


class secondError(RdmlError):
    """Just to have, not used yet"""
    def __init__(self, message):
        RdmlError.__init__(self, message)
    pass


def _get_first_child(base, tag):
    """Get a child element of the base node with a given tag.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements group tag used to select the elements. (string)

    Returns:
        The first child lxml node element found or None.
    """

    for node in base:
        if node.tag.replace("{http://www.rdml.org}", "") == tag:
            return node
    return None


def _get_first_child_text(base, tag):
    """Get a child element of the base node with a given tag.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements group tag used to select the elements. (string)

    Returns:
        The text of first child node element found or an empty string.
    """

    for node in base:
        if node.tag.replace("{http://www.rdml.org}", "") == tag:
            return node.text
    return ""


def _get_first_child_bool(base, tag, triple=True):
    """Get a child element of the base node with a given tag.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements group tag used to select the elements. (string)
        triple: If True, None is returned if not found, if False, False

    Returns:
        The a bool value of tag or if triple is True None.
    """

    for node in base:
        if node.tag.replace("{http://www.rdml.org}", "") == tag:
            return _string_to_bool(node.text, triple)
    if triple is False:
        return False
    else:
        return None


def _string_to_bool(value, triple=True):
    """Translates a string into bool value or None.

    Args:
        value: The string value to evaluate. (string)
        triple: If True, None is returned if not found, if False, False

    Returns:
        The a bool value of tag or if triple is True None.
    """

    if value is None or value == "":
        if triple is True:
            return None
        else:
            return False
    if type(value) is bool:
        return value
    if type(value) is int:
        if value != 0:
            return True
        else:
            return False
    if type(value) is str:
        if value.lower() in ['false', '0', 'f', '-', 'n', 'no']:
            return False
        else:
            return True
    return


def _value_to_booldic(value):
    """Translates a string, list or dic to a dictionary with true/false.

    Args:
        value: The string value to evaluate. (string)

    Returns:
        The a bool value of tag or if triple is True None.
    """

    ret = {}
    if type(value) is str:
        ret[value] = True
    if type(value) is list:
        for ele in value:
            ret[ele] = True
    if type(value) is dict:
        for key, val in value.items():
            ret[key] = _string_to_bool(val, triple=False)
    return ret


def _get_first_child_by_pos_or_id(base, tag, by_id, by_pos):
    """Get a child element of the base node with a given tag and position or id.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements group tag used to select the elements. (string)
        by_id: The unique id to search for. (string)
        by_pos: The position of the element in the list (int)

    Returns:
        The child node element found or raise error.
    """

    if by_id is None and by_pos is None:
        raise RdmlError('Either an ' + tag + ' id or a position must be provided.')
    if by_id is not None and by_pos is not None:
        raise RdmlError('Only an ' + tag + ' id or a position can be provided.')
    allChildren = _get_all_children(base, tag)
    if by_id is not None:
        for node in allChildren:
            if node.get('id') == by_id:
                return node
        raise RdmlError('The ' + tag + ' id: ' + by_id + ' was not found in RDML file.')
    if by_pos is not None:
        if by_pos < 0 or by_pos > len(allChildren) - 1:
            raise RdmlError('The ' + tag + ' position ' + by_pos + ' is out of range.')
        return allChildren[by_pos]


def _add_first_child_to_dic(base, dic, opt, tag):
    """Adds the first child element with a given tag to a dictionary.

    Args:
        base: The base node element. (lxml node)
        dic: The dictionary to add the element to (dictionary)
        opt: If false and id is not found in base, the element is added with an empty string (Bool)
        tag: Child elements group tag used to select the elements. (string)

    Returns:
        The dictionary with the added element.
    """

    for node in base:
        if node.tag.replace("{http://www.rdml.org}", "") == tag:
            dic[tag] = node.text
            return dic
    if not opt:
        dic[tag] = ""
    return dic


def _get_all_children(base, tag):
    """Get a list of all child elements with a given tag.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements group tag used to select the elements. (string)

    Returns:
        A list with all child node elements found or an empty list.
    """

    ret = []
    for node in base:
        if node.tag.replace("{http://www.rdml.org}", "") == tag:
            ret.append(node)
    return ret


def _get_all_children_id(base, tag):
    """Get a list of ids of all child elements with a given tag.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements group tag used to select the elements. (string)

    Returns:
        A list with all child id strings found or an empty list.
    """

    ret = []
    for node in base:
        if node.tag.replace("{http://www.rdml.org}", "") == tag:
            ret.append(node.get('id'))
    return ret


def _get_number_of_children(base, tag):
    """Count all child elements with a given tag.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements group tag used to select the elements. (string)

    Returns:
        A int number of the found child elements with the id.
    """

    counter = 0
    for node in base:
        if node.tag.replace("{http://www.rdml.org}", "") == tag:
            counter += 1
    return counter


def _check_unique_id(base, tag, id):
    """Find all child elements with a given group and check if the id is already used.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements group tag used to select the elements. (string)
        id: The unique id to search for. (string)

    Returns:
        False if the id is already used, True if not.
    """

    for node in base:
        if node.tag.replace("{http://www.rdml.org}", "") == tag:
            if node.get('id') == id:
                return False
    return True


def _create_new_element(base, tag, id):
    """Create a new element with a given tag and id.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements group tag. (string)
        id: The unique id of the new element. (string)

    Returns:
        False if the id is already used, True if not.
    """

    if id is None or id == "":
        raise RdmlError('An ' + tag + ' id must be provided.')
    if not _check_unique_id(base, tag, id):
        raise RdmlError('The ' + tag + ' id "' + id + '" must be unique.')

    return ET.Element(tag, id=id)


def _add_new_subelement(base, basetag, tag, text, opt):
    """Create a new element with a given tag and id.

    Args:
        base: The base node element. (lxml node)
        basetag: Child elements group tag. (string)
        tag: Child elements own tag, to be created. (string)
        text: The text content of the new element. (string)
        opt: If true, the element is optional (Bool)

    Returns:
        Nothing, the base lxml element is modified.
    """

    if opt is False:
        if text is None or text == "":
            raise RdmlError('An ' + basetag + ' ' + tag + ' must be provided.')
        ET.SubElement(base, tag).text = text
    else:
        if text is not None and text != "":
            ET.SubElement(base, tag).text = text


def _change_subelement(base, tag, xmlkeys, value, opt, vtype, id_as_element=False):
    """Change the value of the element with a given tag.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements own tag, to be created. (string)
        xmlkeys: The list of possible keys in the right order for xml (list strings)
        value: The text content of the new element.
        opt: If true, the element is optional (Bool)
        vtype: If true, the element is optional ("string", "int", "float")
        id_as_element: If true, handle tag "id" as element, else as attribute

    Returns:
        Nothing, the base lxml element is modified.
    """

    # Todo validate values with vtype
    goodVal = value
    if vtype == "bool":
        ev = _string_to_bool(value, triple=True)
        if ev is None or ev == "":
            goodVal = ""
        else:
            if ev:
                goodVal = "true"
            else:
                goodVal = "false"

    if opt is False:
        if goodVal is None or goodVal == "":
            raise RdmlError('A value for ' + tag + ' must be provided.')

    if tag == "id" and id_as_element is False:
        if base.get('id') != goodVal:
            par = base.getparent()
            groupTag = base.tag.replace("{http://www.rdml.org}", "")
            if not _check_unique_id(par, groupTag, goodVal):
                raise RdmlError('The ' + groupTag + ' id "' + goodVal + '" is not unique.')
            base.attrib['id'] = goodVal
        return

    # Check if the tag already excists
    elem = _get_first_child(base, tag)
    if elem is not None:
        if goodVal is None or goodVal == "":
            base.remove(elem)
        else:
            elem.text = goodVal
    else:
        if goodVal is not None and goodVal != "":
            new_node = ET.Element(tag)
            new_node.text = goodVal
            place = _get_tag_pos(base, tag, xmlkeys, 0)
            base.insert(place, new_node)


def _get_or_create_subelement(base, tag, xmlkeys):
    """Get element with a given tag, if not present, create it.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements own tag, to be created. (string)
        xmlkeys: The list of possible keys in the right order for xml (list strings)

    Returns:
        The node element with the tag.
    """

    # Check if the tag already excists
    if _get_first_child(base, tag) is None:
        new_node = ET.Element(tag)
        place = _get_tag_pos(base, tag, xmlkeys, 0)
        base.insert(place, new_node)
    return _get_first_child(base, tag)


def _remove_irrelevant_subelement(base, tag):
    """If element with a given tag has no children, remove it.

    Args:
        base: The base node element. (lxml node)
        tag: Child elements own tag, to be created. (string)

    Returns:
        The node element with the tag.
    """

    # Check if the tag already excists
    elem = _get_first_child(base, tag)
    if elem is None:
        return
    if len(elem) == 0:
        base.remove(elem)


def _move_subelement(base, tag, id, xmlkeys, position):
    """Change the value of the element with a given tag.

    Args:
        base: The base node element. (lxml node)
        tag: The id to search for. (string)
        id: The unique id of the new element. (string)
        xmlkeys: The list of possible keys in the right order for xml (list strings)
        position: the new position of the element (int)

    Returns:
        Nothing, the base lxml element is modified.
    """

    pos = _get_tag_pos(base, tag, xmlkeys, position)
    ele = _get_first_child_by_pos_or_id(base, tag, id, None)
    base.insert(pos, ele)


def _get_tag_pos(base, tag, xmlkeys, pos):
    """Returns a position were to add a subelement with the given tag inc. pos offset.

    Args:
        base: The base node element. (lxml node)
        tag: The id to search for. (string)
        xmlkeys: The list of possible keys in the right order for xml (list strings)
        pos: The position relative to the tag elements (int)

    Returns:
        The int number of were to add the element with the tag.
    """

    count = _get_number_of_children(base, tag)
    offset = pos
    if pos is None or pos < 0:
        offset = 0
    if pos > count:
        offset = count
    return _get_first_tag_pos(base, tag, xmlkeys) + offset


def _get_first_tag_pos(base, tag, xmlkeys):
    """Returns a position were to add a subelement with the given tag.

    Args:
        base: The base node element. (lxml node)
        tag: The id to search for. (string)
        xmlkeys: The list of possible keys in the right order for xml (list strings)

    Returns:
        The int number of were to add the element with the tag.
    """

    listrest = xmlkeys[xmlkeys.index(tag):]
    counter = 0
    for node in base:
        if node.tag.replace("{http://www.rdml.org}", "") in listrest:
            return counter
        counter += 1
    return counter


class Rdml:
    """RDML-Python library
    
    The root element used to open, write, read and edit RDML files.
    
    Attributes:
        _rdmlData: The RDML XML object from lxml.
        _node: The root node of the RDML XML object.
        _rdmlVersion: A string like '1.2' with the version of the rdmlData object.
    """

    def __init__(self, filename=None):
        """Inits an empty RDML instance with new() or load RDML file with load().

        Args:
            self: The class self parameter.
            filename: The name of the RDML file to load.

        Returns:
            No return value. Function may raise RdmlError if required.
        """
        
        self._rdmlData = None
        self._node = None
        self._rdmlVersion = '0.0'
        if filename:
            self.load(filename)
        else:
            self.new()

    def __getitem__(self, key):
        """Returns data of the key.

        Args:
            self: The class self parameter.
            key: The key of the experimenter subelement

        Returns:
            A string of the data or None.
        """
        if key == "version":
            return self.version()
        if key in ["dateMade", "dateUpdated"]:
            return _get_first_child_text(self._node, key)
        raise KeyError

    def keys(self):
        """Returns a list of the keys.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["version", "dateMade", "dateUpdated"]

    def xmlkeys(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["dateMade", "dateUpdated", "id", "experimenter", "documentation", "dye",
                "sample", "target", "thermalCyclingConditions", "experiment"]

    def new(self):
        """Creates an new empty RDML object with the current date.

        Args:
            self: The class self parameter.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        data = "<rdml version='1.2' xmlns:rdml='http://www.rdml.org' xmlns='http://www.rdml.org'>\n<dateMade>"
        data += datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
        data += "</dateMade>\n<dateUpdated>"
        data += datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
        data += "</dateUpdated>\n</rdml>"
        self.loadXMLString(data)
        return

    def load(self, filename):
        """Load an RDML file with decompression of rdml_data.xml or an XML file. Uses loadXMLString().

        Args:
            self: The class self parameter.
            filename: The name of the RDML file to load.

        Returns:
            No return value. Function may raise RdmlError if required.
        """    

        if zipfile.is_zipfile(filename):
            zf = zipfile.ZipFile(filename, 'r')
            try:
                data = zf.read('rdml_data.xml').decode('utf-8')
            except KeyError:
                raise RdmlError('No rdml_data.xml in compressed RDML file found.')
            else:
                self.loadXMLString(data)
        else:
            with open(filename, 'r') as txtfile:
                data = txtfile.read()
                if data:
                    self.loadXMLString(data)
                else:
                    raise RdmlError('File format error, not a valid RDML or XML file.')

    def save(self, filename):
        """Save an RDML file with compression of rdml_data.xml.

        Args:
            self: The class self parameter.
            filename: The name of the RDML file to save to.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        elem = _get_first_child(self._node, "dateUpdated")
        elem.text = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
        data = ET.tostring(self._rdmlData, pretty_print=True)
        zf2 = zipfile.ZipFile(filename, mode='w', compression=zipfile.ZIP_DEFLATED,)
        try:
            zf2.writestr('rdml_data.xml', data)
        finally:
            zf2.close()

    def loadXMLString(self, data):
        """Create RDML object from xml string. !ENTITY and DOCSTRINGS will be removed.

        Args:
            self: The class self parameter.
            data: The xml string of the RDML file to load.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        # To avoid some xml attacs based on
        # <!ENTITY entityname "replacement text">
        data = re.sub(r"<\W*!ENTITY[^>]+>", "", data)
        data = re.sub(r"!ENTITY", "", data)
        try:
            self._rdmlData = ET.ElementTree(ET.fromstring(data))
        except ET.XMLSyntaxError:
            raise RdmlError('XML load error, not a valid RDML or XML file.')
        self._node = self._rdmlData.getroot()
        if self._node.tag.replace("{http://www.rdml.org}", "") != 'rdml':
            raise RdmlError('Root element is not \'rdml\', not a valid RDML or XML file.')
        self._rdmlVersion = self._node.get('version')
        # Remainder: Update version in new() and validate()
        if not self._rdmlVersion in ['1.0','1.1','1.2']:
            raise RdmlError('Unknown or unsupported RDML file version.')

    def validate(self, filename=None):
        """Validate the RDML object against its schema or load file and validate it.

        Args:
            self: The class self parameter.
            filename: The name of the RDML file to load.

        Returns:
            A string with the validation result as a two column table.
        """

        notes = ""
        if filename:
            try:
                vd = Rdml(filename)
            except RdmlError as err:
                notes += 'RDML file structure:\tFalse\t' + str(err) + '\n'
                return notes
            notes += "RDML file structure:\tTrue\tValid file structure.\n"
        else:
            vd = self
        version = vd.version()
        rdmlws = os.path.dirname(os.path.abspath(__file__))
        if version == '1.0':
            xmlschema_doc = ET.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_0_REC.xsd'))
        elif version == '1.1':
            xmlschema_doc = ET.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_1_REC.xsd'))
        elif version == '1.2':
            xmlschema_doc = ET.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_2_REC.xsd'))
        else:
            notes += 'RDML version:\tFalse\tUnknown schema version' + version + '\n'
            return notes
        notes += "RDML version:\tTrue\t" + version + "\n"

        xmlschema = ET.XMLSchema(xmlschema_doc)
        result = xmlschema.validate(vd._rdmlData)
        if result:
            notes += 'Schema validation result:\tTrue\tRDML file is valid.\n'
        else:
            notes += 'Schema validation result:\tFalse\tRDML file is not valid.\n'
        log = xmlschema.error_log
        for err in log:
            notes += 'Schema validation error:\tFalse\t'
            notes += "Line %s, Column %s: %s \n" % (err.line, err.column, err.message)
        return notes

    def isvalid(self, filename=None):
        """Validate the RDML object against its schema or load file and validate it.

        Args:
            self: The class self parameter.
            filename: The name of the RDML file to load.

        Returns:
            True or false as the validation result.
        """

        if filename:
            try:
                vd = Rdml(filename)
            except RdmlError as err:
                return False
        else:
            vd = self
        version = vd.version()
        rdmlws = os.path.dirname(os.path.abspath(__file__))
        if version == '1.0':
            xmlschema_doc = ET.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_0_REC.xsd'))
        elif version == '1.1':
            xmlschema_doc = ET.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_1_REC.xsd'))
        elif version == '1.2':
            xmlschema_doc = ET.parse(os.path.join(rdmlws, 'schema', 'RDML_v1_2_REC.xsd'))
        else:
            return False
        xmlschema = ET.XMLSchema(xmlschema_doc)
        result = xmlschema.validate(vd._rdmlData)
        if result:
            return True
        else:
            return False

    def version(self):
        """Returns the version string of the RDML object.

        Args:
            self: The class self parameter.

        Returns:
            A string of the version like '1.1'.
        """

        return self._rdmlVersion

    def rdmlids(self):
        """Returns a list of all rdml id elements.

        Args:
            self: The class self parameter.

        Returns:
            A list of all rdml id elements.
        """

        exp = _get_all_children(self._node, "id")
        ret = []
        for node in exp:
            ret.append(Rdmlid(node, self._rdmlVersion))
        return ret

    def new_rdmlid(self, publisher, serialNumber, MD5Hash=None, newposition=None):
        """Creates a new rdml id element.

        Args:
            self: The class self parameter.
            publisher: Publisher who created the serialNumber (required)
            serialNumber: Serial Number for this file provided by publisher (required)
            MD5Hash: A MD5 hash for this file (optional)

        Returns:
            Nothing, changes self.
        """

        new_node = ET.Element("id")
        _add_new_subelement(new_node, "id", "publisher", publisher, False)
        _add_new_subelement(new_node, "id", "serialNumber", serialNumber, False)
        _add_new_subelement(new_node, "id", "MD5Hash", MD5Hash, True)
        place = _get_tag_pos(self._node, "id", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def move_rdmlid(self, oldposition, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "id", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "id", None, oldposition)
        self._node.insert(pos, ele)

    def get_rdmlid(self, byposition=None):
        """Returns an experimenter element by position or id.

        Args:
            self: The class self parameter.
            byposition: Select the element by position in the list.

        Returns:
            The found element or None.
        """

        return Rdmlid(_get_first_child_by_pos_or_id(self._node, "id", None, byposition),
                      self._rdmlVersion)

    def delete_rdmlid(self, byposition=None):
        """Deletes an experimenter element.

        Args:
            self: The class self parameter.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "id", None, byposition)
        self._node.remove(elem)

    def experimenters(self):
        """Returns a list of all experimenter elements.

        Args:
            self: The class self parameter.

        Returns:
            A list of all experimenter elements.
        """

        exp = _get_all_children(self._node, "experimenter")
        ret = []
        for node in exp:
            ret.append(Experimenter(node, self._rdmlVersion))
        return ret

    def new_experimenter(self, id, firstName, lastName, email=None, labName=None, labAddress=None, newposition=None):
        """Creates a new experimenter element.

        Args:
            self: The class self parameter.
            id: Experimenter unique id
            firstName: Experimenters first name (required)
            lastName: Experimenters last name (required)
            email: Experimenters email (optional)
            labName: Experimenters lab name (optional)
            labAddress: Experimenters lab address (optional)
            newposition: Experimenters position in the list of experimenters (optional)

        Returns:
            Nothing, changes self.
        """
        new_node = _create_new_element(self._node, "experimenter", id)
        _add_new_subelement(new_node, "experimenter", "firstName", firstName, False)
        _add_new_subelement(new_node, "experimenter", "lastName", lastName, False)
        _add_new_subelement(new_node, "experimenter", "email", email, True)
        _add_new_subelement(new_node, "experimenter", "labName", labName, True)
        _add_new_subelement(new_node, "experimenter", "labAddress", labAddress, True)
        place = _get_tag_pos(self._node, "experimenter", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def move_experimenter(self, id, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            id: Experimenter unique id
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        _move_subelement(self._node, "experimenter", id, self.xmlkeys(), newposition)

    def get_experimenter(self, byid=None, byposition=None):
        """Returns an experimenter element by position or id.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            The found element or None.
        """

        return Experimenter(_get_first_child_by_pos_or_id(self._node, "experimenter", byid, byposition),
                            self._rdmlVersion)

    def delete_experimenter(self, byid=None, byposition=None):
        """Deletes an experimenter element.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "experimenter", byid, byposition)
        self._node.remove(elem)
        # Todo delete in all use places

    def documentations(self):
        """Returns a list of all documentation elements.

        Args:
            self: The class self parameter.

        Returns:
            A list of all documentation elements.
        """

        exp = _get_all_children(self._node, "documentation")
        ret = []
        for node in exp:
            ret.append(Documentation(node, self._rdmlVersion))
        return ret

    def new_documentation(self, id, text=None, newposition=None):
        """Creates a new documentation element.

        Args:
            self: The class self parameter.
            id: Documentation unique id
            text: Documentation descriptive test (optional)
            newposition: Experimenters position in the list of experimenters (optional)

        Returns:
            Nothing, changes self.
        """
        new_node = _create_new_element(self._node, "documentation", id)
        _add_new_subelement(new_node, "documentation", "text", text, True)
        place = _get_tag_pos(self._node, "documentation", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def move_documentation(self, id, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            id: Documentation unique id
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        _move_subelement(self._node, "documentation", id, self.xmlkeys(), newposition)

    def get_documentation(self, byid=None, byposition=None):
        """Returns an documentation element by position or id.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            The found element or None.
        """

        return Documentation(_get_first_child_by_pos_or_id(self._node, "documentation", byid, byposition),
                            self._rdmlVersion)

    def delete_documentation(self, byid=None, byposition=None):
        """Deletes an documentation element.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "documentation", byid, byposition)
        self._node.remove(elem)
        # Todo delete in all use places

    def dyes(self):
        """Returns a list of all dye elements.

        Args:
            self: The class self parameter.

        Returns:
            A list of all dye elements.
        """

        exp = _get_all_children(self._node, "dye")
        ret = []
        for node in exp:
            ret.append(Dye(node, self._rdmlVersion))
        return ret

    def new_dye(self, id, description=None, newposition=None):
        """Creates a new dye element.

        Args:
            self: The class self parameter.
            id: Dye unique id
            description: Dye descriptive test (optional)
            newposition: Dye position in the list of dyes (optional)

        Returns:
            Nothing, changes self.
        """
        new_node = _create_new_element(self._node, "dye", id)
        _add_new_subelement(new_node, "dye", "description", description, True)
        place = _get_tag_pos(self._node, "dye", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def move_dye(self, id, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            id: Dye unique id
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        _move_subelement(self._node, "dye", id, self.xmlkeys(), newposition)

    def get_dye(self, byid=None, byposition=None):
        """Returns an dye element by position or id.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            The found element or None.
        """

        return Dye(_get_first_child_by_pos_or_id(self._node, "dye", byid, byposition),
                   self._rdmlVersion)

    def delete_dye(self, byid=None, byposition=None):
        """Deletes an dye element.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "dye", byid, byposition)
        self._node.remove(elem)
        # Todo delete in all use places

    def samples(self):
        """Returns a list of all sample elements.

        Args:
            self: The class self parameter.

        Returns:
            A list of all sample elements.
        """

        exp = _get_all_children(self._node, "sample")
        ret = []
        for node in exp:
            ret.append(Sample(node, self._rdmlVersion))
        return ret

    def new_sample(self, id, type, newposition=None):
        """Creates a new sample element.

        Args:
            self: The class self parameter.
            id: Sample unique id (required)
            type: Sample type (required)
            newposition: Experimenters position in the list of experimenters (optional)

        Returns:
            Nothing, changes self.
        """

        if type not in ["unkn", "ntc", "nac", "std", "ntp", "nrt", "pos", "opt"]:
            raise RdmlError('Unknown or unsupported sample type value "' + type + '".')
        new_node = _create_new_element(self._node, "sample", id)
        _add_new_subelement(new_node, "sample", "type", type, False)
        place = _get_tag_pos(self._node, "sample", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def move_sample(self, id, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            id: Sample unique id
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        _move_subelement(self._node, "sample", id, self.xmlkeys(), newposition)

    def get_sample(self, byid=None, byposition=None):
        """Returns an sample element by position or id.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            The found element or None.
        """

        return Sample(_get_first_child_by_pos_or_id(self._node, "sample", byid, byposition),
                      self._rdmlVersion)

    def delete_sample(self, byid=None, byposition=None):
        """Deletes an sample element.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "sample", byid, byposition)
        self._node.remove(elem)
        # Todo delete in all use places

    def targets(self):
        """Returns a list of all target elements.

        Args:
            self: The class self parameter.

        Returns:
            A list of all target elements.
        """

        exp = _get_all_children(self._node, "target")
        ret = []
        for node in exp:
            ret.append(Target(node, self._rdmlVersion))
        return ret

    def new_target(self, id, type, newposition=None):
        """Creates a new target element.

        Args:
            self: The class self parameter.
            id: Target unique id (required)
            type: Target type (required)
            newposition: Targets position in the list of targets (optional)

        Returns:
            Nothing, changes self.
        """

        if type not in ["ref", "toi"]:
            raise RdmlError('Unknown or unsupported target type value "' + type + '".')
        new_node = _create_new_element(self._node, "target", id)
        _add_new_subelement(new_node, "target", "type", type, False)
        place = _get_tag_pos(self._node, "target", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def move_target(self, id, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            id: Target unique id
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        _move_subelement(self._node, "target", id, self.xmlkeys(), newposition)

    def get_target(self, byid=None, byposition=None):
        """Returns an target element by position or id.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            The found element or None.
        """

        return Target(_get_first_child_by_pos_or_id(self._node, "target", byid, byposition),
                      self._rdmlVersion)

    def delete_target(self, byid=None, byposition=None):
        """Deletes an target element.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "target", byid, byposition)
        self._node.remove(elem)
        # Todo delete in all use places

    def therm_cyc_cons(self):
        """Returns a list of all thermalCyclingConditions elements.

        Args:
            self: The class self parameter.

        Returns:
            A list of all target elements.
        """

        exp = _get_all_children(self._node, "thermalCyclingConditions")
        ret = []
        for node in exp:
            ret.append(Therm_cyc_cons(node, self._rdmlVersion))
        return ret

    def new_therm_cyc_cons(self, id, newposition=None):
        """Creates a new thermalCyclingConditions element.

        Args:
            self: The class self parameter.
            id: ThermalCyclingConditions unique id (required)
            newposition: Targets position in the list of targets (optional)

        Returns:
            Nothing, changes self.
        """

        new_node = _create_new_element(self._node, "thermalCyclingConditions", id)
        step = ET.SubElement(new_node, "step")
        ET.SubElement(step, "nr").text = "1"
        ET.SubElement(step, "lidOpen")
        place = _get_tag_pos(self._node, "thermalCyclingConditions", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def move_therm_cyc_cons(self, id, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            id: ThermalCyclingConditions unique id
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        _move_subelement(self._node, "thermalCyclingConditions", id, self.xmlkeys(), newposition)

    def get_therm_cyc_cons(self, byid=None, byposition=None):
        """Returns an thermalCyclingConditions element by position or id.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            The found element or None.
        """

        return Therm_cyc_cons(_get_first_child_by_pos_or_id(self._node, "thermalCyclingConditions", byid, byposition),
                              self._rdmlVersion)

    def delete_therm_cyc_cons(self, byid=None, byposition=None):
        """Deletes an thermalCyclingConditions element.

        Args:
            self: The class self parameter.
            byid: Select the element by the element id.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "thermalCyclingConditions", byid, byposition)
        self._node.remove(elem)
        # Todo delete in all use places

    def tojson(self):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """

        allRdmlids = self.rdmlids()
        rdmlids = []
        for elem in allRdmlids:
            rdmlids.append(elem.tojson())

        allExperimenters = self.experimenters()
        experimenters = []
        for exp in allExperimenters:
            experimenters.append(exp.tojson())

        allDocumentations = self.documentations()
        documentations = []
        for exp in allDocumentations:
            documentations.append(exp.tojson())

        allDyes = self.dyes()
        dyes = []
        for exp in allDyes:
            dyes.append(exp.tojson())

        allSamples = self.samples()
        samples = []
        for exp in allSamples:
            samples.append(exp.tojson())

        allTargets = self.targets()
        targets = []
        for exp in allTargets:
            targets.append(exp.tojson())

        allTherm_cyc_cons = self.therm_cyc_cons()
        therm_cyc_cons = []
        for exp in allTherm_cyc_cons:
            therm_cyc_cons.append(exp.tojson())

        data = {
            "rdml": {
                "version": self["version"],
                "dateMade": self["dateMade"],
                "dateUpdated": self["dateUpdated"],
                "ids": rdmlids,
                "experimenters": experimenters,
                "documentations": documentations,
                "dyes": dyes,
                "samples": samples,
                "targets": targets,
                "therm_cyc_cons": therm_cyc_cons,
                "experiments": []
            }
        }
        return data


class Rdmlid:
    """RDML-Python library

    The rdml id element used to read and edit one experimenter.

    Attributes:
        _node: The rdml id node of the RDML XML object.
        _rdmlVersion: A string like '1.2' with the version of the rdmlData object.
    """

    def __init__(self, node, version):
        """Inits an rdml id instance.

        Args:
            self: The class self parameter.
            node: The experimenter node.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        self._node = node
        self._rdmlVersion = version

    def __getitem__(self, key):
        """Returns the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the experimenter subelement

        Returns:
            A string of the data or None.
        """

        if key in ["publisher", "serialNumber"]:
            return _get_first_child_text(self._node, key)
        if key in ["MD5Hash"]:
            var = _get_first_child_text(self._node, key)
            if var == "":
                return None
            else:
                return var
        raise KeyError

    def __setitem__(self, key, value):
        """Changes the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the experimenter subelement
            value: The new value for the key

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """
        if key in ["publisher", "serialNumber"]:
            return _change_subelement(self._node, key, self.xmlkeys(), value, False, "string")
        if key in ["MD5Hash"]:
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "string")
        raise KeyError

    def keys(self):
        """Returns a list of the keys.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["publisher", "serialNumber", "MD5Hash"]

    def xmlkeys(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return self.keys()

    def tojson(self):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """

        data = {
            "publisher": _get_first_child_text(self._node, "publisher"),
            "serialNumber": _get_first_child_text(self._node, "serialNumber")
        }
        _add_first_child_to_dic(self._node, data, True, "MD5Hash")
        return data


class Experimenter:
    """RDML-Python library

    The experimenter element used to read and edit one experimenter.

    Attributes:
        _node: The experimenter node of the RDML XML object.
        _rdmlVersion: A string like '1.2' with the version of the rdmlData object.
    """

    def __init__(self, node, version):
        """Inits an experimenter instance.

        Args:
            self: The class self parameter.
            node: The experimenter node.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        self._node = node
        self._rdmlVersion = version

    def __getitem__(self, key):
        """Returns the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the experimenter subelement

        Returns:
            A string of the data or None.
        """

        if key == "id":
            return self._node.get('id')
        if key in ["firstName", "lastName"]:
            return _get_first_child_text(self._node, key)
        if key in ["email", "labName", "labAddress"]:
            var = _get_first_child_text(self._node, key)
            if var == "":
                return None
            else:
                return var
        raise KeyError

    def __setitem__(self, key, value):
        """Changes the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the experimenter subelement
            value: The new value for the key

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """
        if key in ["id", "firstName", "lastName"]:
            return _change_subelement(self._node, key, self.xmlkeys(), value, False, "string")
        if key in ["email", "labName", "labAddress"]:
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "string")
        raise KeyError

    def keys(self):
        """Returns a list of the keys.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["id", "firstName", "lastName", "email", "labName", "labAddress"]

    def xmlkeys(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return self.keys()

    def tojson(self):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """

        data = {
            "id": self._node.get('id'),
            "firstName": _get_first_child_text(self._node, "firstName"),
            "lastName": _get_first_child_text(self._node, "lastName")
        }
        _add_first_child_to_dic(self._node, data, True, "email")
        _add_first_child_to_dic(self._node, data, True, "labName")
        _add_first_child_to_dic(self._node, data, True, "labAddress")
        return data


class Documentation:
    """RDML-Python library

    The documentation element used to read and edit one documentation tag.

    Attributes:
        _node: The documentation node of the RDML XML object.
        _rdmlVersion: A string like '1.2' with the version of the rdmlData object.
    """

    def __init__(self, node, version):
        """Inits an documentation instance.

        Args:
            self: The class self parameter.
            node: The documentation node.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        self._node = node
        self._rdmlVersion = version

    def __getitem__(self, key):
        """Returns the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the documentation subelement

        Returns:
            A string of the data or None.
        """

        if key == "id":
            return self._node.get('id')
        if key == "text":
            var = _get_first_child_text(self._node, key)
            if var == "":
                return None
            else:
                return var
        raise KeyError

    def __setitem__(self, key, value):
        """Changes the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the documentation subelement
            value: The new value for the key

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """
        if key == "id":
            return _change_subelement(self._node, key, self.xmlkeys(), value, False, "string")
        if key == "text":
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "string")
        raise KeyError

    def keys(self):
        """Returns a list of the keys.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["id", "text"]

    def xmlkeys(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return self.keys()

    def tojson(self):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """

        data = {
            "id": self._node.get('id'),
        }
        _add_first_child_to_dic(self._node, data, True, "text")
        return data


class Dye:
    """RDML-Python library

    The dye element used to read and edit one dye.

    Attributes:
        _node: The dye node of the RDML XML object.
        _rdmlVersion: A string like '1.2' with the version of the rdmlData object.
    """

    def __init__(self, node, version):
        """Inits an dye instance.

        Args:
            self: The class self parameter.
            node: The dye node.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        self._node = node
        self._rdmlVersion = version

    def __getitem__(self, key):
        """Returns the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the dye subelement

        Returns:
            A string of the data or None.
        """

        if key == "id":
            return self._node.get('id')
        if key == "description":
            var = _get_first_child_text(self._node, key)
            if var == "":
                return None
            else:
                return var
        raise KeyError

    def __setitem__(self, key, value):
        """Changes the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the dye subelement
            value: The new value for the key

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """
        if key == "id":
            return _change_subelement(self._node, key, self.xmlkeys(), value, False, "string")
        if key == "description":
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "string")
        raise KeyError

    def keys(self):
        """Returns a list of the keys.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["id", "description"]

    def xmlkeys(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return self.keys()

    def tojson(self):
        """Returns a json of the RDML object.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """

        data = {
            "id": self._node.get('id'),
        }
        _add_first_child_to_dic(self._node, data, True, "description")
        return data


class Sample:
    """RDML-Python library

    The samples element used to read and edit one sample.

    Attributes:
        _node: The sample node of the RDML XML object.
        _rdmlVersion: A string like '1.2' with the version of the rdmlData object.
    """

    def __init__(self, node, version):
        """Inits an sample instance.

        Args:
            self: The class self parameter.
            node: The sample node.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        self._node = node
        self._rdmlVersion = version

    def __getitem__(self, key):
        """Returns the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the sample subelement

        Returns:
            A string of the data or None.
        """

        if key == "id":
            return self._node.get('id')
        if key == "type":
            return _get_first_child_text(self._node, key)
        if key == "description":
            var = _get_first_child_text(self._node, key)
            if var == "":
                return None
            else:
                return var
        if key in ["interRunCalibrator", "calibratorSample"]:
            return _get_first_child_bool(self._node, key, triple=True)
        if key in ["quantity", "templateRNAQuantity", "templateDNAQuantity"]:
            ele = _get_first_child(self._node, key)
            vdic = {}
            vdic["value"] = _get_first_child_text(ele, "value")
            vdic["unit"] = _get_first_child_text(ele, "unit")
            if len(vdic.keys()) != 0:
                return vdic
            else:
                return None
        if key in ["templateRNAQuality", "templateDNAQuality"]:
            ele = _get_first_child(self._node, key)
            vdic = {}
            vdic["method"] = _get_first_child_text(ele, "method")
            vdic["result"] = _get_first_child_text(ele, "result")
            if len(vdic.keys()) != 0:
                return vdic
            else:
                return None
        if key in ["cdnaSynthesisMethod_enzyme", "cdnaSynthesisMethod_primingMethod",
                   "cdnaSynthesisMethod_dnaseTreatment", "cdnaSynthesisMethod_thermalCyclingConditions"]:
            ele = _get_first_child(self._node, "cdnaSynthesisMethod")
            if ele is None:
                return None
            if key == "cdnaSynthesisMethod_enzyme":
                return _get_first_child_text(ele, "enzyme")
            if key == "cdnaSynthesisMethod_primingMethod":
                return _get_first_child_text(ele, "primingMethod")
            if key == "cdnaSynthesisMethod_dnaseTreatment":
                return _get_first_child_text(ele, "dnaseTreatment")
            if key == "cdnaSynthesisMethod_thermalCyclingConditions":
                forId = _get_first_child(ele, "thermalCyclingConditions")
                if forId is not None:
                    return forId.attrib['id']
                else:
                    return None
            raise RdmlError('Sample cdnaSynthesisMethod programming read error.')
        raise KeyError

    def __setitem__(self, key, value):
        """Changes the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the sample subelement
            value: The new value for the key

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        if key == "type":
            if value not in ["unkn", "ntc", "nac", "std", "ntp", "nrt", "pos", "opt"]:
                raise RdmlError('Unknown or unsupported sample type value "' + value + '".')

        if key in ["id", "type"]:
            return _change_subelement(self._node, key, self.xmlkeys(), value, False, "string")
        if key == "description":
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "string")
        if key in ["interRunCalibrator", "calibratorSample"]:
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "bool")
        if key in ["quantity", "templateRNAQuantity", "templateDNAQuantity"]:
            if value is None:
                return
            if "value" not in value or "unit" not in value:
                raise RdmlError('Sample ' + key + ' must have a dictionary with "value" and "unit" as value.')
            if value["unit"] not in ["", "cop", "fold", "dil", "ng", "nMol", "other"]:
                raise RdmlError('Unknown or unsupported sample ' + key + ' value "' + value + '".')
            ele = _get_or_create_subelement(self._node, key, self.xmlkeys())
            _change_subelement(ele, "value", ["value", "unit"], value["value"], True, "float")
            if value["value"] != "":
                _change_subelement(ele, "unit", ["value", "unit"], value["unit"], True, "string")
            else:
                _change_subelement(ele, "unit", ["value", "unit"], "", True, "string")
            _remove_irrelevant_subelement(self._node, key)
            return
        if key in ["templateRNAQuality", "templateDNAQuality"]:
            if value is None:
                return
            if "method" not in value or "result" not in value:
                raise RdmlError('"' + key + '" must have a dictionary with "method" and "result" as value.')
            ele = _get_or_create_subelement(self._node, key, self.xmlkeys())
            _change_subelement(ele, "method", ["method", "result"], value["method"], True, "string")
            _change_subelement(ele, "result", ["method", "result"], value["result"], True, "float")
            _remove_irrelevant_subelement(self._node, key)
            return
        if key in ["cdnaSynthesisMethod_enzyme", "cdnaSynthesisMethod_primingMethod",
                   "cdnaSynthesisMethod_dnaseTreatment", "cdnaSynthesisMethod_thermalCyclingConditions"]:
            ele = _get_or_create_subelement(self._node, "cdnaSynthesisMethod", self.xmlkeys())
            if key == "cdnaSynthesisMethod_enzyme":
                _change_subelement(ele, "enzyme",
                                   ["enzyme", "primingMethod", "dnaseTreatment", "thermalCyclingConditions"],
                                   value, True, "string")
            if key == "cdnaSynthesisMethod_primingMethod":
                if value not in ["", "oligo-dt", "random", "target-specific", "oligo-dt and random", "other"]:
                    raise RdmlError('Unknown or unsupported sample ' + key + ' value "' + value + '".')
                _change_subelement(ele, "primingMethod",
                                   ["enzyme", "primingMethod", "dnaseTreatment", "thermalCyclingConditions"],
                                   value, True, "string")
            if key == "cdnaSynthesisMethod_dnaseTreatment":
                _change_subelement(ele, "dnaseTreatment",
                                   ["enzyme", "primingMethod", "dnaseTreatment", "thermalCyclingConditions"],
                                   value, True, "bool")
            if key == "cdnaSynthesisMethod_thermalCyclingConditions":
                forId = _get_or_create_subelement(ele, "thermalCyclingConditions",
                                                  ["enzyme", "primingMethod", "dnaseTreatment",
                                                   "thermalCyclingConditions"])
                if value is not None and value != "":
                    # Todo check ID
                    forId.attrib['id'] = value
                else:
                    ele.remove(forId)
            _remove_irrelevant_subelement(self._node, "cdnaSynthesisMethod")
            return
        raise KeyError

    def keys(self):
        """Returns a list of the keys.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["id", "description", "type", "interRunCalibrator", "quantity", "calibratorSample",
                "cdnaSynthesisMethod_enzyme", "cdnaSynthesisMethod_primingMethod",
                "cdnaSynthesisMethod_dnaseTreatment", "cdnaSynthesisMethod_thermalCyclingConditions",
                "templateRNAQuantity", "templateRNAQuality", "templateDNAQuantity", "templateDNAQuality"]

    def xmlkeys(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["description", "documentation", "xRef", "type", "interRunCalibrator",
                "quantity", "calibratorSample", "cdnaSynthesisMethod",
                "templateRNAQuantity", "templateRNAQuality", "templateDNAQuantity", "templateDNAQuality"]

    def xrefs(self):
        """Returns a list of the xrefs in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of dics with name and id strings.
        """

        xref = _get_all_children(self._node, "xRef")
        ret = []
        for node in xref:
            data = {}
            _add_first_child_to_dic(node, data, True, "name")
            _add_first_child_to_dic(node, data, True, "id")
            ret.append(data)
        return ret

    def new_xref(self, name=None, id=None, newposition=None):
        """Creates a new xrefs element.

        Args:
            self: The class self parameter.
            name: Publisher who created the xRef
            id: Serial Number for this sample provided by publisher
            newposition: The new position of the element

        Returns:
            Nothing, changes self.
        """

        if name is None and id is None:
            raise RdmlError('Either name or id is required to create a xRef.')
        new_node = ET.Element("xRef")
        _add_new_subelement(new_node, "xRef", "name", name, True)
        _add_new_subelement(new_node, "xRef", "id", id, True)
        place = _get_tag_pos(self._node, "xRef", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def edit_xref(self, oldposition, newposition=None, name=None, id=None):
        """Creates a new xrefs element.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element
            name: Publisher who created the xRef
            id: Serial Number for this sample provided by publisher

        Returns:
            Nothing, changes self.
        """

        if oldposition is None:
            raise RdmlError('A oldposition is required to edit a xRef.')
        if (name is None or name == "") and (id is None or id == ""):
            self.delete_xref(oldposition)
            return
        pos = _get_tag_pos(self._node, "xRef", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "xRef", None, oldposition)
        _change_subelement(ele, "name", ["name", "id"], name, True, "string")
        _change_subelement(ele, "id", ["name", "id"], id, True, "string", id_as_element=True)
        self._node.insert(pos, ele)

    def move_xref(self, oldposition, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "xRef", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "xRef", None, oldposition)
        self._node.insert(pos, ele)

    def delete_xref(self, byposition):
        """Deletes an experimenter element.

        Args:
            self: The class self parameter.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "xRef", None, byposition)
        self._node.remove(elem)

    def documentation_ids(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return _get_all_children_id(self._node, "documentation")

    def update_documentation_ids(self, ids):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.
            ids: A dictionary with id and true/false pairs

        Returns:
            True if a change was made, else false. Function may raise RdmlError if required.
        """

        old = self.documentation_ids()
        good_ids = _value_to_booldic(ids)
        mod = False

        for id, inc in good_ids.items():
            if inc is True:
                if id not in old:
                    new_node = _create_new_element(self._node, "documentation", id)
                    place = _get_tag_pos(self._node, "documentation", self.xmlkeys(), 999999999)
                    self._node.insert(place, new_node)
                    mod = True
            else:
                if id in old:
                    elem = _get_first_child_by_pos_or_id(self._node, "documentation", id, None)
                    self._node.remove(elem)
                    mod = True
        return mod

    def move_documentation(self, oldposition, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "documentation", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "documentation", None, oldposition)
        self._node.insert(pos, ele)

    def tojson(self):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """

        data = {
            "id": self._node.get('id'),
        }
        _add_first_child_to_dic(self._node, data, True, "description")
        data["documentations"] = self.documentation_ids()
        data["xRefs"] = self.xrefs()
        _add_first_child_to_dic(self._node, data, False, "type")
        _add_first_child_to_dic(self._node, data, True, "interRunCalibrator")
        elem = _get_first_child(self._node, "quantity")
        if elem is not None:
            qdic = {}
            _add_first_child_to_dic(elem, qdic, False, "value")
            _add_first_child_to_dic(elem, qdic, False, "unit")
            data["quantity"] = qdic
        _add_first_child_to_dic(self._node, data, True, "calibratorSample")
        elem = _get_first_child(self._node, "cdnaSynthesisMethod")
        if elem is not None:
            qdic = {}
            _add_first_child_to_dic(elem, qdic, True, "enzyme")
            _add_first_child_to_dic(elem, qdic, True, "primingMethod")
            _add_first_child_to_dic(elem, qdic, True, "dnaseTreatment")
            forId = _get_first_child(elem, "thermalCyclingConditions")
            if forId is not None:
                if forId.attrib['id'] != "":
                    qdic["thermalCyclingConditions"] = forId.attrib['id']
            if len(qdic.keys()) != 0:
                data["cdnaSynthesisMethod"] = qdic
        elem = _get_first_child(self._node, "templateRNAQuantity")
        if elem is not None:
            qdic = {}
            _add_first_child_to_dic(elem, qdic, False, "value")
            _add_first_child_to_dic(elem, qdic, False, "unit")
            data["templateRNAQuantity"] = qdic
        elem = _get_first_child(self._node, "templateRNAQuality")
        if elem is not None:
            qdic = {}
            _add_first_child_to_dic(elem, qdic, False, "method")
            _add_first_child_to_dic(elem, qdic, False, "result")
            data["templateRNAQuality"] = qdic
        elem = _get_first_child(self._node, "templateDNAQuantity")
        if elem is not None:
            qdic = {}
            _add_first_child_to_dic(elem, qdic, False, "value")
            _add_first_child_to_dic(elem, qdic, False, "unit")
            data["templateDNAQuantity"] = qdic
        elem = _get_first_child(self._node, "templateDNAQuality")
        if elem is not None:
            qdic = {}
            _add_first_child_to_dic(elem, qdic, False, "method")
            _add_first_child_to_dic(elem, qdic, False, "result")
            data["templateDNAQuality"] = qdic
        return data


class Target:
    """RDML-Python library

    The target element used to read and edit one target.

    Attributes:
        _node: The target node of the RDML XML object.
        _rdmlVersion: A string like '1.2' with the version of the rdmlData object.
    """

    def __init__(self, node, version):
        """Inits an target instance.

        Args:
            self: The class self parameter.
            node: The target node.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        self._node = node
        self._rdmlVersion = version

    def __getitem__(self, key):
        """Returns the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the target subelement

        Returns:
            A string of the data or None.
        """

        if key == "id":
            return self._node.get('id')
        if key == "type":
            return _get_first_child_text(self._node, key)
        if key in ["description", "amplificationEfficiencyMethod", "amplificationEfficiency", "detectionLimit"]:
            var = _get_first_child_text(self._node, key)
            if var == "":
                return None
            else:
                return var
        if key == "dyeId":
            forId = _get_first_child(self._node, key)
            if forId is not None:
                return forId.attrib['id']
            else:
                return None
        if key in ["sequences_forwardPrimer_threePrimeTag", "sequences_forwardPrimer_fivePrimeTag",
                   "sequences_forwardPrimer_sequence", "sequences_reversePrimer_threePrimeTag",
                   "sequences_reversePrimer_fivePrimeTag", "sequences_reversePrimer_sequence",
                   "sequences_probe1_threePrimeTag", "sequences_probe1_fivePrimeTag",
                   "sequences_probe1_sequence", "sequences_probe2_threePrimeTag",
                   "sequences_probe2_fivePrimeTag", "sequences_probe2_sequence",
                   "sequences_amplicon_threePrimeTag", "sequences_amplicon_fivePrimeTag",
                   "sequences_amplicon_sequence"]:
            prim = _get_first_child(self._node, "sequences")
            if prim is None:
                return None
            sec = None
            if key in ["sequences_forwardPrimer_threePrimeTag", "sequences_forwardPrimer_fivePrimeTag",
                       "sequences_forwardPrimer_sequence"]:
                sec = _get_first_child(prim, "forwardPrimer")
            if key in ["sequences_reversePrimer_threePrimeTag", "sequences_reversePrimer_fivePrimeTag",
                       "sequences_reversePrimer_sequence"]:
                sec = _get_first_child(prim, "reversePrimer")
            if key in ["sequences_probe1_threePrimeTag", "sequences_probe1_fivePrimeTag", "sequences_probe1_sequence"]:
                sec = _get_first_child(prim, "probe1")
            if key in ["sequences_probe2_threePrimeTag", "sequences_probe2_fivePrimeTag", "sequences_probe2_sequence"]:
                sec = _get_first_child(prim, "probe2")
            if key in ["sequences_amplicon_threePrimeTag", "sequences_amplicon_fivePrimeTag",
                       "sequences_amplicon_sequence"]:
                sec = _get_first_child(prim, "amplicon")
            if sec is None:
                return None
            if key in ["sequences_forwardPrimer_threePrimeTag", "sequences_reversePrimer_threePrimeTag",
                       "sequences_probe1_threePrimeTag", "sequences_probe2_threePrimeTag",
                       "sequences_amplicon_threePrimeTag"]:
                return _get_first_child_text(sec, "threePrimeTag")
            if key in ["sequences_forwardPrimer_fivePrimeTag", "sequences_reversePrimer_fivePrimeTag",
                       "sequences_probe1_fivePrimeTag", "sequences_probe2_fivePrimeTag",
                       "sequences_amplicon_fivePrimeTag"]:
                return _get_first_child_text(sec, "fivePrimeTag")
            if key in ["sequences_forwardPrimer_sequence", "sequences_reversePrimer_sequence",
                       "sequences_probe1_sequence", "sequences_probe2_sequence",
                       "sequences_amplicon_sequence"]:
                return _get_first_child_text(sec, "sequence")
            raise RdmlError('Target sequences programming read error.')
        if key in ["commercialAssay_company", "commercialAssay_orderNumber"]:
            prim = _get_first_child(self._node, "commercialAssay")
            if prim is None:
                return None
            if key == "commercialAssay_company":
                return _get_first_child_text(prim, "company")
            if key == "commercialAssay_orderNumber":
                return _get_first_child_text(prim, "orderNumber")
        raise KeyError

    def __setitem__(self, key, value):
        """Changes the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the target subelement
            value: The new value for the key

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        if key == "type":
            if value not in ["ref", "toi"]:
                raise RdmlError('Unknown or unsupported target type value "' + value + '".')

        if key in ["id", "type"]:
            return _change_subelement(self._node, key, self.xmlkeys(), value, False, "string")
        if key in ["description", "amplificationEfficiencyMethod"]:
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "string")
        if key in ["amplificationEfficiency", "detectionLimit"]:
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "float")
        if key == "dyeId":
            forId = _get_or_create_subelement(self._node, "dyeId", self.xmlkeys())
            if value is not None and value != "":
                # Todo check ID
                forId.attrib['id'] = value
            else:
                self._node.remove(forId)
            return
        if key in ["sequences_forwardPrimer_threePrimeTag", "sequences_forwardPrimer_fivePrimeTag",
                   "sequences_forwardPrimer_sequence", "sequences_reversePrimer_threePrimeTag",
                   "sequences_reversePrimer_fivePrimeTag", "sequences_reversePrimer_sequence",
                   "sequences_probe1_threePrimeTag", "sequences_probe1_fivePrimeTag",
                   "sequences_probe1_sequence", "sequences_probe2_threePrimeTag",
                   "sequences_probe2_fivePrimeTag", "sequences_probe2_sequence",
                   "sequences_amplicon_threePrimeTag", "sequences_amplicon_fivePrimeTag",
                   "sequences_amplicon_sequence"]:
            prim = _get_or_create_subelement(self._node, "sequences", self.xmlkeys())
            sec = None
            if key in ["sequences_forwardPrimer_threePrimeTag", "sequences_forwardPrimer_fivePrimeTag",
                       "sequences_forwardPrimer_sequence"]:
                sec = _get_or_create_subelement(prim, "forwardPrimer",
                                                ["forwardPrimer", "reversePrimer", "probe1", "probe2", "amplicon"])
            if key in ["sequences_reversePrimer_threePrimeTag", "sequences_reversePrimer_fivePrimeTag",
                       "sequences_reversePrimer_sequence"]:
                sec = _get_or_create_subelement(prim, "reversePrimer",
                                                ["forwardPrimer", "reversePrimer", "probe1", "probe2", "amplicon"])
            if key in ["sequences_probe1_threePrimeTag", "sequences_probe1_fivePrimeTag", "sequences_probe1_sequence"]:
                sec = _get_or_create_subelement(prim, "probe1",
                                                ["forwardPrimer", "reversePrimer", "probe1", "probe2", "amplicon"])
            if key in ["sequences_probe2_threePrimeTag", "sequences_probe2_fivePrimeTag", "sequences_probe2_sequence"]:
                sec = _get_or_create_subelement(prim, "probe2",
                                                ["forwardPrimer", "reversePrimer", "probe1", "probe2", "amplicon"])
            if key in ["sequences_amplicon_threePrimeTag", "sequences_amplicon_fivePrimeTag",
                       "sequences_amplicon_sequence"]:
                sec = _get_or_create_subelement(prim, "amplicon",
                                                ["forwardPrimer", "reversePrimer", "probe1", "probe2", "amplicon"])
            if sec is None:
                return None
            if key in ["sequences_forwardPrimer_threePrimeTag", "sequences_reversePrimer_threePrimeTag",
                       "sequences_probe1_threePrimeTag", "sequences_probe2_threePrimeTag",
                       "sequences_amplicon_threePrimeTag"]:
                _change_subelement(sec, "threePrimeTag",
                                   ["threePrimeTag", "fivePrimeTag", "sequence"], value, True, "string")
            if key in ["sequences_forwardPrimer_fivePrimeTag", "sequences_reversePrimer_fivePrimeTag",
                       "sequences_probe1_fivePrimeTag", "sequences_probe2_fivePrimeTag",
                       "sequences_amplicon_fivePrimeTag"]:
                _change_subelement(sec, "fivePrimeTag",
                                   ["threePrimeTag", "fivePrimeTag", "sequence"], value, True, "string")
            if key in ["sequences_forwardPrimer_sequence", "sequences_reversePrimer_sequence",
                       "sequences_probe1_sequence", "sequences_probe2_sequence",
                       "sequences_amplicon_sequence"]:
                _change_subelement(sec, "sequence",
                                   ["threePrimeTag", "fivePrimeTag", "sequence"], value, True, "string")
            if key in ["sequences_forwardPrimer_threePrimeTag", "sequences_forwardPrimer_fivePrimeTag",
                       "sequences_forwardPrimer_sequence"]:
                _remove_irrelevant_subelement(prim, "forwardPrimer")
            if key in ["sequences_reversePrimer_threePrimeTag", "sequences_reversePrimer_fivePrimeTag",
                       "sequences_reversePrimer_sequence"]:
                _remove_irrelevant_subelement(prim, "reversePrimer")
            if key in ["sequences_probe1_threePrimeTag", "sequences_probe1_fivePrimeTag", "sequences_probe1_sequence"]:
                _remove_irrelevant_subelement(prim, "probe1")
            if key in ["sequences_probe2_threePrimeTag", "sequences_probe2_fivePrimeTag", "sequences_probe2_sequence"]:
                _remove_irrelevant_subelement(prim, "probe2")
            if key in ["sequences_amplicon_threePrimeTag", "sequences_amplicon_fivePrimeTag",
                       "sequences_amplicon_sequence"]:
                _remove_irrelevant_subelement(prim, "amplicon")
            _remove_irrelevant_subelement(self._node, "sequences")
            return
        if key in ["commercialAssay_company", "commercialAssay_orderNumber"]:
            ele = _get_or_create_subelement(self._node, "commercialAssay", self.xmlkeys())
            if key == "commercialAssay_company":
                _change_subelement(ele, "company", ["company", "orderNumber"], value, True, "string")
            if key == "commercialAssay_orderNumber":
                _change_subelement(ele, "orderNumber", ["company", "orderNumber"], value, True, "string")
            _remove_irrelevant_subelement(self._node, "commercialAssay")
            return
        raise KeyError

    def keys(self):
        """Returns a list of the keys.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["id", "description", "type", "amplificationEfficiencyMethod", "amplificationEfficiency",
                "detectionLimit", "dyeId", "sequences_forwardPrimer_threePrimeTag",
                "sequences_forwardPrimer_fivePrimeTag", "sequences_forwardPrimer_sequence",
                "sequences_reversePrimer_threePrimeTag", "sequences_reversePrimer_fivePrimeTag",
                "sequences_reversePrimer_sequence", "sequences_probe1_threePrimeTag",
                "sequences_probe1_fivePrimeTag", "sequences_probe1_sequence", "sequences_probe2_threePrimeTag",
                "sequences_probe2_fivePrimeTag", "sequences_probe2_sequence", "sequences_amplicon_threePrimeTag",
                "sequences_amplicon_fivePrimeTag", "sequences_amplicon_sequence", "commercialAssay_company",
                "commercialAssay_orderNumber"]

    def xmlkeys(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["description", "documentation", "xRef", "type", "amplificationEfficiencyMethod",
                "amplificationEfficiency", "detectionLimit", "dyeId", "sequences", "commercialAssay"]

    def xrefs(self):
        """Returns a list of the xrefs in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of dics with name and id strings.
        """

        xref = _get_all_children(self._node, "xRef")
        ret = []
        for node in xref:
            data = {}
            _add_first_child_to_dic(node, data, True, "name")
            _add_first_child_to_dic(node, data, True, "id")
            ret.append(data)
        return ret

    def new_xref(self, name=None, id=None, newposition=None):
        """Creates a new xrefs element.

        Args:
            self: The class self parameter.
            name: Publisher who created the xRef
            id: Serial Number for this target provided by publisher
            newposition: The new position of the element

        Returns:
            Nothing, changes self.
        """

        if name is None and id is None:
            raise RdmlError('Either name or id is required to create a xRef.')
        new_node = ET.Element("xRef")
        _add_new_subelement(new_node, "xRef", "name", name, True)
        _add_new_subelement(new_node, "xRef", "id", id, True)
        place = _get_tag_pos(self._node, "xRef", self.xmlkeys(), newposition)
        self._node.insert(place, new_node)

    def edit_xref(self, oldposition, newposition=None, name=None, id=None):
        """Creates a new xrefs element.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element
            name: Publisher who created the xRef
            id: Serial Number for this target provided by publisher

        Returns:
            Nothing, changes self.
        """

        if oldposition is None:
            raise RdmlError('A oldposition is required to edit a xRef.')
        if (name is None or name == "") and (id is None or id == ""):
            self.delete_xref(oldposition)
            return
        pos = _get_tag_pos(self._node, "xRef", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "xRef", None, oldposition)
        _change_subelement(ele, "name", ["name", "id"], name, True, "string")
        _change_subelement(ele, "id", ["name", "id"], id, True, "string", id_as_element=True)
        self._node.insert(pos, ele)

    def move_xref(self, oldposition, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "xRef", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "xRef", None, oldposition)
        self._node.insert(pos, ele)

    def delete_xref(self, byposition):
        """Deletes an experimenter element.

        Args:
            self: The class self parameter.
            byposition: Select the element by position in the list.

        Returns:
            Nothing, changes self.
        """

        elem = _get_first_child_by_pos_or_id(self._node, "xRef", None, byposition)
        self._node.remove(elem)

    def documentation_ids(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return _get_all_children_id(self._node, "documentation")

    def update_documentation_ids(self, ids):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.
            ids: A dictionary with id and true/false pairs

        Returns:
            True if a change was made, else false. Function may raise RdmlError if required.
        """

        old = self.documentation_ids()
        good_ids = _value_to_booldic(ids)
        mod = False

        for id, inc in good_ids.items():
            if inc is True:
                if id not in old:
                    new_node = _create_new_element(self._node, "documentation", id)
                    place = _get_tag_pos(self._node, "documentation", self.xmlkeys(), 999999999)
                    self._node.insert(place, new_node)
                    mod = True
            else:
                if id in old:
                    elem = _get_first_child_by_pos_or_id(self._node, "documentation", id, None)
                    self._node.remove(elem)
                    mod = True
        return mod

    def move_documentation(self, oldposition, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "documentation", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "documentation", None, oldposition)
        self._node.insert(pos, ele)

    def tojson(self):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """

        data = {
            "id": self._node.get('id'),
        }
        _add_first_child_to_dic(self._node, data, True, "description")
        data["documentations"] = self.documentation_ids()
        data["xRefs"] = self.xrefs()
        _add_first_child_to_dic(self._node, data, False, "type")
        _add_first_child_to_dic(self._node, data, True, "amplificationEfficiencyMethod")
        _add_first_child_to_dic(self._node, data, True, "amplificationEfficiency")
        _add_first_child_to_dic(self._node, data, True, "detectionLimit")
        forId = _get_first_child(self._node, "dyeId")
        if forId is not None:
            if forId.attrib['id'] != "":
                data["dyeId"] = forId.attrib['id']
        elem = _get_first_child(self._node, "sequences")
        if elem is not None:
            qdic = {}
            sec = _get_first_child(elem, "forwardPrimer")
            if sec is not None:
                sdic = {}
                _add_first_child_to_dic(sec, sdic, True, "threePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "fivePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "sequence")
                if len(sdic.keys()) != 0:
                    qdic["forwardPrimer"] = sdic
            sec = _get_first_child(elem, "reversePrimer")
            if sec is not None:
                sdic = {}
                _add_first_child_to_dic(sec, sdic, True, "threePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "fivePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "sequence")
                if len(sdic.keys()) != 0:
                    qdic["reversePrimer"] = sdic
            sec = _get_first_child(elem, "probe1")
            if sec is not None:
                sdic = {}
                _add_first_child_to_dic(sec, sdic, True, "threePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "fivePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "sequence")
                if len(sdic.keys()) != 0:
                    qdic["probe1"] = sdic
            sec = _get_first_child(elem, "probe2")
            if sec is not None:
                sdic = {}
                _add_first_child_to_dic(sec, sdic, True, "threePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "fivePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "sequence")
                if len(sdic.keys()) != 0:
                    qdic["probe2"] = sdic
            sec = _get_first_child(elem, "amplicon")
            if sec is not None:
                sdic = {}
                _add_first_child_to_dic(sec, sdic, True, "threePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "fivePrimeTag")
                _add_first_child_to_dic(sec, sdic, True, "sequence")
                if len(sdic.keys()) != 0:
                    qdic["amplicon"] = sdic
            if len(qdic.keys()) != 0:
                data["sequences"] = qdic
        elem = _get_first_child(self._node, "commercialAssay")
        if elem is not None:
            qdic = {}
            _add_first_child_to_dic(elem, qdic, True, "company")
            _add_first_child_to_dic(elem, qdic, True, "orderNumber")
            if len(qdic.keys()) != 0:
                data["commercialAssay"] = qdic
        return data


class Therm_cyc_cons:
    """RDML-Python library

    The thermalCyclingConditions element used to read and edit one thermal Cycling Conditions.

    Attributes:
        _node: The thermalCyclingConditions node of the RDML XML object.
        _rdmlVersion: A string like '1.2' with the version of the rdmlData object.
    """

    def __init__(self, node, version):
        """Inits an thermalCyclingConditions instance.

        Args:
            self: The class self parameter.
            node: The thermalCyclingConditions node.

        Returns:
            No return value. Function may raise RdmlError if required.
        """

        self._node = node
        self._rdmlVersion = version

    def __getitem__(self, key):
        """Returns the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the thermalCyclingConditions subelement

        Returns:
            A string of the data or None.
        """

        if key == "id":
            return self._node.get('id')
        if key in ["description", "lidTemperature"]:
            var = _get_first_child_text(self._node, key)
            if var == "":
                return None
            else:
                return var
        raise KeyError

    def __setitem__(self, key, value):
        """Changes the value for the key.

        Args:
            self: The class self parameter.
            key: The key of the thermalCyclingConditions subelement
            value: The new value for the key

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        if key == "id":
            return _change_subelement(self._node, key, self.xmlkeys(), value, False, "string")
        if key == "description":
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "string")
        if key == "lidTemperature":
            return _change_subelement(self._node, key, self.xmlkeys(), value, True, "float")
        raise KeyError

    def keys(self):
        """Returns a list of the keys.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["id", "description", "lidTemperature"]

    def xmlkeys(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return ["description", "documentation", "lidTemperature", "experimenter", "step"]

    def documentation_ids(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return _get_all_children_id(self._node, "documentation")

    def update_documentation_ids(self, ids):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.
            ids: A dictionary with id and true/false pairs

        Returns:
            True if a change was made, else false. Function may raise RdmlError if required.
        """

        old = self.documentation_ids()
        good_ids = _value_to_booldic(ids)
        mod = False

        for id, inc in good_ids.items():
            if inc is True:
                if id not in old:
                    new_node = _create_new_element(self._node, "documentation", id)
                    place = _get_tag_pos(self._node, "documentation", self.xmlkeys(), 999999999)
                    self._node.insert(place, new_node)
                    mod = True
            else:
                if id in old:
                    elem = _get_first_child_by_pos_or_id(self._node, "documentation", id, None)
                    self._node.remove(elem)
                    mod = True
        return mod

    def move_documentation(self, oldposition, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "documentation", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "documentation", None, oldposition)
        self._node.insert(pos, ele)

    def experimenter_ids(self):
        """Returns a list of the keys in the xml file.

        Args:
            self: The class self parameter.

        Returns:
            A list of the key strings.
        """

        return _get_all_children_id(self._node, "experimenter")

    def update_experimenter_ids(self, ids):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.
            ids: A dictionary with id and true/false pairs

        Returns:
            True if a change was made, else false. Function may raise RdmlError if required.
        """

        old = self.experimenter_ids()
        good_ids = _value_to_booldic(ids)
        mod = False

        for id, inc in good_ids.items():
            if inc is True:
                if id not in old:
                    new_node = _create_new_element(self._node, "experimenter", id)
                    place = _get_tag_pos(self._node, "experimenter", self.xmlkeys(), 999999999)
                    self._node.insert(place, new_node)
                    mod = True
            else:
                if id in old:
                    elem = _get_first_child_by_pos_or_id(self._node, "experimenter", id, None)
                    self._node.remove(elem)
                    mod = True
        return mod

    def move_experimenter(self, oldposition, newposition):
        """Moves the element to the new position in the list.

        Args:
            self: The class self parameter.
            oldposition: The old position of the element
            newposition: The new position of the element

        Returns:
            No return value, changes self. Function may raise RdmlError if required.
        """

        pos = _get_tag_pos(self._node, "experimenter", self.xmlkeys(), newposition)
        ele = _get_first_child_by_pos_or_id(self._node, "experimenter", None, oldposition)
        self._node.insert(pos, ele)

    def tojson(self):
        """Returns a json of the RDML object without fluorescence data.

        Args:
            self: The class self parameter.

        Returns:
            A json of the data.
        """

        data = {
            "id": self._node.get('id'),
        }
        _add_first_child_to_dic(self._node, data, True, "description")
        data["documentations"] = self.documentation_ids()
        _add_first_child_to_dic(self._node, data, True, "lidTemperature")
        data["experimenters"] = self.experimenter_ids()
        return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='The command line interface to the RDML-Python library.')
    parser.add_argument('-v', '--validate', metavar="data.rdml", dest='validate', help='validate file against schema')
    parser.add_argument("-d", "--doooo", dest="doooo", help="just do stuff")

    args = parser.parse_args()

    # Validate RDML file
    if args.validate:
        inst = Rdml()
        res = inst.validate(filename=args.validate)
        print(res)
        sys.exit(0)

    # Tryout things
    if args.doooo:
        print('Tryout')
        xx = Rdml('rdml_data.xml')
        xx.save('new.rdml')
