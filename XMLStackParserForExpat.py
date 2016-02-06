# spud/spud/utils/XMLStackParserForExpat.py
# A parser and Expat XML parser helper object that allows decorator-based path event handlers.
import xml.parsers.expat
import re
from collections import defaultdict
import threading

# module global regex objects...
COMPONENT_REGEX = re.compile(r"/?[^/]+")
ATTRIB_VALUE_REGEX = re.compile("\[\s*\@([a-zA-z0-9]+)\s*\=\s*\'(.*)\'\s*\]")
ATTRIB_NOVALUE_REGEX = re.compile("\[\s*\@([a-zA-z0-9]+)\s*\]")
ATTRIB_CAPTURE_REGEX = re.compile("\[\s*\@([a-zA-z0-9]+)\s*\=\s*\$([a-zA-z0-9]+)\s*\]")
ATTRIB_TEXT_REGEX = re.compile("\[\s*\^text\s*\=\s*\$([a-zA-z0-9]+)\s*\]")

class XMLExpatParserHelperPredicate(object):
    def evaluate(self, attrs_stack):
        """
            Needs to be implemented in derived classes.
            Should return a tuple of (status, attribute_captures, text_captures)
            where the status is the True/False result of testing the predicate,
            and attribute_captures and text_captures are None or each are a
            dict that contains the name -> value mapping for the captured
            attribute or text values. See the path expression capture patterns
            documented in XMLExpatParserHelper and XMLStackParser.
        """
        raise NotImplementedError
    
class XMLExpatParserHelperAttributeValuePredicate(XMLExpatParserHelperPredicate):
    def __init__(self, index, attrib, value):
        self.index = index
        self.attrib = attrib
        self.value = value
    
    def evaluate(self, attrs_stack, text_stack):
        index = self.index
        attrib = self.attrib
        value = self.value
        if attrib not in attrs_stack[index]:
            return (False, None, None)
        if value != attrs_stack[index][attrib]:
            return (False, None, None)
        return (True, None, None)
    
class XMLExpatParserHelperAttributePresentPredicate(XMLExpatParserHelperPredicate):
    def __init__(self, index, attrib):
        self.index = index
        self.attrib = attrib

    def evaluate(self, attrs_stack, text_stack):
        if self.attrib not in attrs_stack[self.index]:
            return (False, None, None)
        return (True, None, None)

class XMLExpatParserHelperAttributeCapturePredicate(XMLExpatParserHelperPredicate):
    def __init__(self, index, attrib, symbol):
        self.index = index
        self.attrib = attrib
        self.symbol = symbol
    
    def evaluate(self, attrs_stack, text_stack):
        index = self.index
        attrib = self.attrib
        if attrib not in attrs_stack[index]:
            return (False, None, None)
        return (True, {self.symbol : attrs_stack[index][attrib]}, None)

class XMLExpatParserHelperTextCapturePredicate(XMLExpatParserHelperPredicate):
    def __init__(self, index, symbol):
        self.index = index
        self.symbol = symbol
    
    def evaluate(self, attrs_stack, text_stack):
        index = self.index
        return (True, None, {self.symbol : text_stack[index]})

def parse_xpath(xpath):
    """
        Parse a subset of the xpath expression syntax into a plainpath part and a list of predicates part.
        plainpath is the xpath without any predicate expressions.
        predicates is a list of predicate tuples of the form (attribute, value, index), where index is the
        stack offset of the element in the predicate path to which the predicate condition attribute == value
        applies.
        
        An index of -1 indicates the top of the stack.
        A value of None indicates that the predicate tests for the attribute presence, but not its value.
        Return (plainpath, predicates)
    """
    index = 0
    components = []
    predicates = []
    terms = COMPONENT_REGEX.findall(xpath)
    nterms = len(terms)
    for term in terms:
        index += 1
        component = ATTRIB_VALUE_REGEX.sub("", term)
        component = ATTRIB_NOVALUE_REGEX.sub("", component)
        component = ATTRIB_CAPTURE_REGEX.sub("", component)
        component = ATTRIB_TEXT_REGEX.sub("", component)
        components.append(component)
        for predicate in ATTRIB_VALUE_REGEX.findall(term):
            (attrib, value) = predicate
            predicates.append(XMLExpatParserHelperAttributeValuePredicate(index - nterms - 1, attrib, value))
        for predicate in ATTRIB_NOVALUE_REGEX.findall(term):
            attrib = predicate
            predicates.append(XMLExpatParserHelperAttributePresentPredicate(index - nterms - 1, attrib))
        for predicate in ATTRIB_CAPTURE_REGEX.findall(term):
            (attrib, symbol) = predicate
            predicates.append(XMLExpatParserHelperAttributeCapturePredicate(index - nterms - 1, attrib, symbol))
        for predicate in ATTRIB_TEXT_REGEX.findall(term):
            symbol = predicate
            predicates.append(XMLExpatParserHelperTextCapturePredicate(index - nterms - 1, symbol))

    plainpath = "".join(components)
    return (plainpath, predicates)

# Callback object wraps methods and binds to instance derived from XMLExpatParserHelper
# This is a place to temporarily store the parsed xpath's after application of the @xml_start_handler
# @xml_end_handler decorators before the XMLExpatParserHelper instance exists.
# The XMLExpatParserHelperElementCallback class type is also a marker for the XMLExpatParserHelper
# instance to know which methods have been wrapped and contain path callback handler information.
class XMLExpatParserHelperElementCallback(object):
    def __init__(self, func):
        # save original function and initialize path lists...
        self.func = func
        self.name = func.__name__
        self.startElementPathPredicatesList = []
        self.endElementPathPredicatesList = []
        
    def __call__(self, *args):
        # call original function with arguments,
        # func is a data attribute, so self is not passed to the call automatically...
        self.func(*args)
        
    def bind(self, target):
        # attach function to target to make callable as method of target...
        def _bind(*args):
            self(target, *args)
        return _bind
    
def xml_start_handler(elementpath):
    """
        Decorates a method of a class derived from XMLExpatParserHelper.
    """
    def _xml_start_handler(func):
        if not isinstance(func, XMLExpatParserHelperElementCallback):
            func = XMLExpatParserHelperElementCallback(func)
        pathpredicates = parse_xpath(elementpath)   
        func.startElementPathPredicatesList.append(pathpredicates)
        return func
    return _xml_start_handler

def xml_end_handler(elementpath):
    def _xml_end_handler(func):
        if not isinstance(func, XMLExpatParserHelperElementCallback):
            func = XMLExpatParserHelperElementCallback(func)
        pathpredicates = parse_xpath(elementpath)   
        func.endElementPathPredicatesList.append(pathpredicates)
        return func
    return _xml_end_handler
    
class AttributeAccessibleDict(dict):
    """
        Same as a regular Python dict, except that existing elements
        of the dictionary can be accessed with the member dot
        notation, instead of the square bracket quoted string notation.
        This will make writing XMLExpatParserHelper and XMLStackParser
        element callback handlers more concise.

        For example dx["foo"] can be accessed as dx.foo, which is clearer.
    """
    def __getattr__(self, name):
        return self[name]
        
class XMLExpatParserHelper(object):
    """
        Helper class to make using XML Expat parser easier.
        
        1. Derive subclass from XMLExpatParserHelper, and add path event methods. The __init__
        constructor should always pass a valid Expat parser object and a user-defined context object.
       
        2. Decorate path event methods with @xml_start_handler and @xml_end_handler.
            Path expressions may be absolute or relative, and may include attribute predicates:
            a. "/ROOT/TOPIC/ARTICLE" - absolute path
            b. "TOPIC/ARTICLE" - relative path
            c. "TOPIC[@name]/ARTICLE" - relative path requiring 'name' attribute on the TOPIC element
            d. "TOPIC[@name='foo']/ARTICLE" - relative path requiring a value of "foo" on the 'name' attribute
            on the TOPIC element
            e. "TOPIC[@name=$foo]/ARTICLE" - relative path requiring a 'name' attribute on the TOPIC element,
            the value of 'name' will be passed to the event method handler as an entry in the captures
            dictionary parameter.
            f. "TOPIC[^text=$foo]/ARTICLE" - relative path that captures the text on the TOPIC element when
            the path is matched and returns all of the text that has been seen so far within the scope of
            the TOPIC element as the 'foo' entry in the texts capture dictionary. Note that ^text is a special
            symbol referring to the text within the predicate's element.
            
            Multiple captures can be present on a path or an individual element in a path.
            Note that only single quotes are supported in path predicates, and the '/' cannot occur in a predicate.         
       
        3. Path event methods must be of the form:
            arbitrary_method_name(self, element_name, element_attributes, captures_dictionary, texts_dictionary, context)
            where: context is the context object passed to the constructor.
            also: the dictionaries element_attributes, captures_dictionary, and texts_dictionary are instances
            of AttributeAccessibleDict, which allows dictionary items to be accessed as attributes using
            the dot notation.
            
        4. Path methods may make use the accumulate_text() method to turn on or off the saving of text elements.     
       
        5. Invoke the parser.Parse() or parser.ParseFile() methods as usual.
        The appropriate methods will be invoked when the elements specified by the
        path expressions with fulfilled predicates start and end.
    """
    def __init__(self, parser, context):
        # set up initial values...
        self.context = context
        parser.StartElementHandler = self.start_element
        parser.EndElementHandler = self.end_element
        parser.CharacterDataHandler = self.char_data
        self.stack = []
        self.attrs_stack = []
        self.text_stack = []        
        self.accumulate = True
        self.startElementPathPredicatesTable = defaultdict(lambda : [])
        self.endElementPathPredicatesTable = defaultdict(lambda : [])
        self.startElementPathPredicatesTableGet = self.startElementPathPredicatesTable.get
        self.endElementPathPredicatesTableGet = self.endElementPathPredicatesTable.get
        # attach decorated class methods to the parser helper,
        # loop over all methods and construct path tables for this object by looking for callable
        # objects that have "start_element_path" and "end_element_path" attributes...
        for o in self.__class__.__dict__.values():
            if isinstance(o, XMLExpatParserHelperElementCallback):
                self.bind_callback_object(o)
    
    def bind_callback_object(self, callback):
        # bind XMLExpatParserHelperElementCallback object __call__ method to the XMLExpatParserHelper...
        bound = callback.bind(self)
        # insert into the appropriate path tables...
        for (path, predicates) in callback.startElementPathPredicatesList:
            self.startElementPathPredicatesTable[path].append((bound, predicates))
        for (path, predicates) in callback.endElementPathPredicatesList:
            self.endElementPathPredicatesTable[path].append((bound, predicates))
        
    # these methods are part of the internal implementation
    # of the XMLExpatParserHelper base class...
    def start_element(self, name, attrs):
        # convert attrs to AttributeAccessibleDict so that callback dictionaries
        # are all of the same type, this requires an extra dict creation for
        # attrs, but not for captures or texts...
        attrs = AttributeAccessibleDict(attrs) 
        self.stack.append(name)
        self.attrs_stack.append(attrs)
        self.text_stack.append([])
        (handler, captures, texts) = self.current_stack_to_start_element_handler()
        handler(name, attrs, captures, texts, self.context)
    
    def end_element(self, name):
        attrs = self.attrs_stack[-1]
        (handler, captures, texts) = self.current_stack_to_end_element_handler()
        handler(name, attrs, captures, texts, self.context)
        self.attrs_stack.pop()
        self.stack.pop()
        subtext = [t for t in self.text_stack.pop() if t]
        if self.text_stack:
            self.text_stack[-1].extend(subtext)
        
    def char_data(self, data):
        if self.accumulate:
            self.text_stack[-1].append(data)
    
    def stack_to_relative_path(self, stack):
        return "/".join(stack)
        
    def stack_to_absolute_path(self, stack):
        return "/%s" % self.stack_to_relative_path(stack)
        
    def get_matching_handler(self, stack, path, pathTableGet):
        # get the handler that matches the path and the predicates...
        handler = None
        captures = None
        texts = None
        for (bound, predicates) in pathTableGet(path, []):
            # find first bound handler for the path where there are no predicates
            # or all predicates are True...
            if not predicates:
                handler = bound
                break
            (status, captures, texts) = self.evaluate_predicates(predicates)
            if status:
                handler = bound
                break
        return (handler, captures, texts)
        
    def current_stack_to_element_handler(self, pathTableGet, no_match_element_handler):
        # copy current stack so that we can modify it...
        stack = self.stack[:]
        # check for absolute path match...
        path = self.stack_to_absolute_path(stack)
        (handler, captures, texts) = self.get_matching_handler(stack, path, pathTableGet)
        # check for relative path subpath matches...
        while stack and not handler:
            path = self.stack_to_relative_path(stack)
            (handler, captures, texts) = self.get_matching_handler(stack, path, pathTableGet)
            stack = stack[1:]
        # clear out values if no handler matched...
        if not handler:
            handler = no_match_element_handler
            captures = None
            texts = None
        elif texts:
            # handler matched, collapse texts lists to text strings in an
            # AttributeAccessibleDict at last possible step...
            texts = AttributeAccessibleDict([(symbol, "".join(txts)) for (symbol, txts) in texts.items()])
        return (handler, captures, texts)

    def evaluate_predicates(self, predicates):
        # check for any False predicates...
        texts = {}
        captures = AttributeAccessibleDict()
        for p in predicates:
            (status, capture, text) = p.evaluate(self.attrs_stack, self.text_stack)
            if not status:
                return (False, None, None)
            else:
                if capture:
                    # accumulate attribute, symbol captures if any...
                    captures.update(capture)
                if text:
                    texts.update(text)
        # all predicates were True...
        return (True, captures, texts)
    
    def current_stack_to_start_element_handler(self):
        return self.current_stack_to_element_handler(self.startElementPathPredicatesTableGet, self.no_match_start_element_handler)
    
    def current_stack_to_end_element_handler(self):
        return self.current_stack_to_element_handler(self.endElementPathPredicatesTableGet, self.no_match_end_element_handler)
    
    # these are default methods called when no handler matches an element start or end.
    # these can be overridden if needed, the base class implementation does nothing...
    def no_match_start_element_handler(self, name, attrs, captures, texts, context):
        return
    
    def no_match_end_element_handler(self, name, attrs, captures, texts, context):
        return

    def accumulate_text(self, onOrOff):
        """
            Turn text accumulation on or off.
            Returns the prior accumulate state.
        """
        old = self.accumulate
        self.accumulate = onOrOff
        return old

        
class XMLStackParser(object):
    """
        XMLStackParser uses the Expat XML parser and the XMLExpatParserHelper to construct
        a simple Python 'domain specific language' for parsing XML files. Usage is similar
        to XMLExpatParserHelper in that function decorators are used to declare the path
        expressions that match a handler. 
        
        With XMLStackParser, the @xml_stack_start_handler, and @xml_stack_end_handler
        decorators are used instead of the ones used with XMLExpatParserHelper.

        Path expressions may be absolute or relative, and may include attribute predicates:
        a. "/ROOT/TOPIC/ARTICLE" - absolute path
        b. "TOPIC/ARTICLE" - relative path
        c. "TOPIC[@name]/ARTICLE" - relative path requiring 'name' attribute on the TOPIC element
        d. "TOPIC[@name='foo']/ARTICLE" - relative path requiring a value of "foo" on the 'name' attribute
        on the TOPIC element
        e. "TOPIC[@name=$foo]/ARTICLE" - relative path requiring a 'name' attribute on the TOPIC element,
        the value of 'name' will be passed to the event method handler as an entry in the captures
        dictionary parameter.
        f. "TOPIC[^text=$foo]/ARTICLE" - relative path that captures the text on the TOPIC element when
        the path is matched and returns all of the text that has been seen so far within the scope of
        the TOPIC element as the 'foo' entry in the texts capture dictionary. Note that ^text is a special
        symbol referring to the text within the predicate's element.
        
        Multiple captures can be present on a path or an individual element in a path.
        Note that only single quotes are supported in path predicates, and the '/' cannot occur in a predicate.         
                
        Path event methods must be of the form:
            arbitrary_method_name(self, element_name, element_attributes, captures_dictionary, texts_dictionary, context)
            where: context is the context object passed to the constructor.
            also: the dictionaries element_attributes, captures_dictionary, and texts_dictionary are instances
            of AttributeAccessibleDict, which allows dictionary items to be accessed as attributes using
            the dot notation.
            
        Path methods may make use the accumulate_text() method to turn on or off the saving of text elements.     
       
        Invoke the parser.parse() method at the end of the with block to perform the parsing.
        The appropriate methods will be invoked when the elements specified by the
        path expressions with fulfilled predicates start and end.
                                
        The XMLStackParser is intended to be used as a context manager object, taking
        advantage of the Python 'with' keyword. In general, the use of XMLStackParser
        will look like this:
        
        context = [] # context should be a mutable data type
        with XMLStackParser(file, context) as parser:
            # declare handler functions and decorate with @xml_stack_start_handler,
            # and @xml_stack_end_handler, for example:
            
            @xml_stack_start_handler(/SOME/LONG/PATH[^text=$txt1]/WITH[^text=$txt2]/PREDICATES[@foo='bar'])
            some_element_handler(helper, name, attrs, captures, texts, context):
                # do something when the path expression is matched, for example grab the texts...
                context.append(" ".join(texts.items())
        
            # after the handlers have been declared and decorated, parse the file:
            parser.parse()
        # at this point the parse is complete, and the context object contains
        # whatever data put into it by the handlers...
        print context
        
        XMLStackParser provides the following methods that can be called within
        the with block:
            parse()
            accumulate_text(yesNo)
            
        Don't forget to invoke the XMLStackParser.parser() method as the last
        step in the with block!
    """
    # class thread local storage...
    thread_local_storage = threading.local()
    
    def __init__(self, textOrFile, context):
        # create Expat parser and save state...
        self.expat = xml.parsers.expat.ParserCreate()
        self.textOrFile = textOrFile
        self.context = context
        self.helper = None
        self.callbackHandlers = []
        
    def __enter__(self):
        # if no parser stack on this thread then add one...
        if not hasattr(XMLStackParser.thread_local_storage, "xmlStackParserStack"):
            XMLStackParser.thread_local_storage.xmlStackParserStack = []
        # set ourselves to be the current XMLStackParser object for this thread...
        XMLStackParser.thread_local_storage.xmlStackParserStack.append(self)        
        # create the parser helper object and attach to the parser...
        self.helper = XMLExpatParserHelper(self.expat, self.context)
        # handler methods will now be attached to the helper object via running the decorators...
        return self
        
    def __exit__(self, type, value, traceback):
        # clear the current XMLStackParser object for this thread...
        XMLStackParser.thread_local_storage.xmlStackParserStack.pop()        
        return False

    @classmethod
    def get_thread_current_stack_parser(cls):
        return XMLStackParser.thread_local_storage.xmlStackParserStack[-1]
    
    def parse(self):
        # attach the handlers to the helper object on the XMLStackParser on the top of the stack...
        for handler in self.callbackHandlers:
            self.get_thread_current_stack_parser().helper.bind_callback_object(handler)  
        # run the parser in the correct mode for the data...
        if isinstance(self.textOrFile, str) or isinstance(self.textOrFile, unicode):
            self.expat.Parse(self.textOrFile)
        else:
            self.expat.ParseFile(self.textOrFile)
            
    def add_path_handler(self, elementpath, callback, isstarthandler):
        """
            Wraps callback with a XMLExpatParserHelperElementCallback
            object if necessary, stores the parsed elementpath in the callback
            object and then attaches the callback to the self
            XMLStackParser object.
            
            Returns a XMLExpatParserHelperElementCallback object that can be used
            to invoke the original function.            
        """
        if not isinstance(callback, XMLExpatParserHelperElementCallback):
            # wrap function as callback object...
            callback = XMLExpatParserHelperElementCallback(callback)
            # add only once to the handler list...
            self.callbackHandlers.append(callback)
        pathpredicates = parse_xpath(elementpath)
        if isstarthandler:
            callback.startElementPathPredicatesList.append(pathpredicates)
        else:
            callback.endElementPathPredicatesList.append(pathpredicates)
        return callback

    def add_path_start_handler(self, elementpath, callback):
        """
            Invokes add_path_handler() with isstarthandler == True.
        
            Returns a XMLExpatParserHelperElementCallback object that can be used
            to invoke the original function.        
        """
        return self.add_path_handler(elementpath, callback, True)

    def add_path_end_handler(self, elementpath, callback):
        """
            Invokes add_path_handler() with isstarthandler == False.
        
            Returns a XMLExpatParserHelperElementCallback object that can be used
            to invoke the original function.        
        """
        return self.add_path_handler(elementpath, callback, False)
                         
    def accumulate_text(self, onOrOff):
        """
            Turn text accumulation on or off.
            Returns the prior accumulate state.
        """
        return self.helper.accumulate_text(onOrOff)

                         
def xml_stack_start_handler(elementpath):
    """
        Decorates a function to use as an element start handler with the current
        XMLStackParser for the active thread.
        
        Calls XMLStackParser.add_path_start_handler().

        Returns a XMLExpatParserHelperElementCallback object that can be used
        to invoke the original function.        
    """
    def _xml_stack_start_handler(func):
        return XMLStackParser.get_thread_current_stack_parser().add_path_start_handler(elementpath, func)
    return _xml_stack_start_handler
    
def xml_stack_end_handler(elementpath):
    """
        Decorates a function to use as an element end handler with the current
        XMLStackParser for the active thread.

        Calls XMLStackParser.add_path_end_handler().
        
        Returns a XMLExpatParserHelperElementCallback object that can be used
        to invoke the original function.
    """
    def _xml_stack_end_handler(func):
        return XMLStackParser.get_thread_current_stack_parser().add_path_end_handler(elementpath, func)
    return _xml_stack_end_handler


