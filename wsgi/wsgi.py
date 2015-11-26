import re, os, time, sys, json, cgi
from cgi import escape
from hashlib import md5
import shutil
from collections import namedtuple
import logging

#from pyvirtualdisplay import Display

WEB_APP_BASE_PATH = "/usr/local/www/CoreTracker/"

sys.path.insert(0,'/anaconda/lib/python2.7/site-packages')

sys.path.insert(0, WEB_APP_BASE_PATH)
sys.path.append(os.path.dirname(__file__))

from Bio.Alphabet import IUPAC
from Bio import AlignIO
from Bio import Alphabet
from Bio.Align import AlignInfo
from server import WebTreeApplication  # Required to use the webplugin
from coretracker import testing_a_lot, coretracker
import utils
from ete2 import PhyloTree
from ete2 import TreeStyle
from ete2 import faces
import posixpath
import threading
from Queue import Queue

TREE = "phylotree.nw"
TMP = WEB_APP_BASE_PATH+"/tmp/"
global IC_CONTENT
IC_CONTENT = []
AALIST = "ACDEFGHIKLMNPQRSTVWY"
global BadRunList
BadRunList = Queue()

ALIGNMENT = {"filtered": "alignment.fasta_filtered", "global" : "alignment.fasta", "ungapped": "alignment.fasta_ungapped", "ID_filtered" : "alignment.fasta_matchremoved"}
EXPIRE_TIME = 60
CHECK_TIME = 10
#tree style 
ts = TreeStyle()
ts.branch_vertical_margin = 5
ts.show_scale = False
THREAD_MAX = 2
#ts.scale = 20

# In order to extend the default WebTreeApplication, we define our own
# WSGI function that handles URL queries
q = Queue()

def call_coretracker(args, setting):
    return testing_a_lot(args, setting)

def del_run(runid=None):
    if runid:
        try:
            folder = os.path.join(WEB_APP_BASE_PATH, 'input', runid)
            shutil.rmtree(folder)
            return True
        except:
            return False
    return False


class ThreadParameters():
    RES = {}
    ERRORS = []
    EXECUTING = {}

ThreadParameters.lock = threading.Lock()

global ThreadParam
ThreadParam = ThreadParameters()


class RunInfo(object):
    """
    This is a class that we keep information for each coretracker 
    run instance
    """
    def __init__(self, name, ctime, result):
        self._id = name
        self._ctime = ctime
        self._result = [str(x) for x in result]

    def get_ctime(self):
        return self._ctime

    def get_result(self):
        return self._result

    def get_id(self):
        return self._id


class CoreTrackerThread(threading.Thread):
    """
    The concept of this class was copied from 
    http://softwareramblings.com/2008/06/running-functions-as-threads-in-python.html
    """

    def __init__(self, target, queue, **kwargs):
        self._target = target
        self.queue = queue
        threading.Thread.__init__(self, **kwargs)

    def run(self):
        global ThreadParam
        global CHECK_TIME
        while True:
            args, setting, runid = self.queue.get(True)
            with ThreadParameters.lock:
                ThreadParameters.EXECUTING[runid] = True
            logging.debug('Starting coretracker execution for runid : '+ runid)
            # hopping this will execute this item
            try :
                answer = self._target(args, setting)
                with ThreadParameters.lock:
                    ThreadParam.RES[runid] = RunInfo(runid, time.time(), answer)
                    del ThreadParam.EXECUTING[runid]

                logging.debug('SUCCESS : Finishing coretracker execution for runid = '+ runid)
            except:
                with ThreadParameters.lock:
                    ThreadParam.ERRORS.append(runid)
                logging.debug('ERROR: Finishing coretracker execution for runid = '+ runid)
            self.queue.task_done()
            # this won't hold the CPU
            time.sleep(CHECK_TIME)

class RunDelThread(threading.Thread):

    def __init__(self, target, wt, ut, queue,  **kwargs):
        self._target = target
        self._wt = wt
        self._ut = ut
        self._queue = queue
        threading.Thread.__init__(self, **kwargs)
        
    def run(self):
        global ThreadParam
        while True:
            time.sleep(self._ut)
            to_delete = []
            with ThreadParameters.lock:
                while len(ThreadParam.ERRORS):
                    to_delete.append(ThreadParameters.ERRORS.pop())
                c_time = time.time()
                for r in ThreadParam.RES.keys():
                    r_info = ThreadParam.RES[r]
                    if c_time - r_info.get_ctime() > self._wt:
                        to_delete.append(r)
                        del ThreadParam.RES[r]
            
            while not self._queue.empty():
                to_delete.append(self._queue.get())

            for runid in to_delete:
                ntry = 3
                deleted = self._target(runid)
                if deleted:
                    logging.debug('Delete run ' + runid)

                else :
                    while not deleted and ntry>0:
                        logging.debug('Unable to delete run ' + runid + ' will try again in 5s')
                        time.sleep(5)
                        deleted = self._target(runid)
                        ntry -= 1
                    if deleted:
                        logging.debug('Was finally able to delete run ' + runid)
            # get in sleep mode and retry in self._ut time
   

def coretracker_run(environ, start_response):
    fieldlist = ['outdir', 'seq', 'tree', 'dnaseq', 'scale', 'excludegap', 'iccontent',\
                'idfilter', 'debug', 'verbose' ]
    CRTArgs = namedtuple("CRTArgs", fieldlist)
    start_response('200 OK', [('Content-Type', 'application/json')])
    res = {}
    if environ['REQUEST_METHOD'].upper() == 'GET' and  environ['QUERY_STRING']:
        queries = cgi.parse_qs(environ['QUERY_STRING'])
        runid = queries.get('time', '')[0]
        folder = os.path.join(WEB_APP_BASE_PATH, 'input', runid)
        if os.path.exists(folder):
            
            args =  CRTArgs(TMP, queries.get('sseq', [None])[0], 
                        queries.get('stree', [None])[0], 
                        queries.get('gseq', [None])[0], 1.0, 
                        float(queries.get('gapfilter')[0]),
                        float(queries.get('icfilter')[0]),
                        float(queries.get('idfilter')[0]) , True, 2
            )
            
            setting = utils.Settings()
            EXCLUDE_AA = "".join([x for x in AALIST if x not in queries.get('includeAA')[0]])
            AA_MAJORITY_THRESH = float(queries.get('aamajthresh')[0])
            FREQUENCY_THRESHOLD = float(queries.get('freqthresh')[0])
            GENETIC_CODE = int(queries.get('gcode')[0])
            COUNT_THRESHOLD = int(queries.get('ctthresh')[0])
            LIMIT_TO_SUSPECTED_SPECIES = queries.get('limsuspected')[0] == 'on'

            setting.set(EXCLUDE_AA=EXCLUDE_AA, AA_MAJORITY_THRESH=AA_MAJORITY_THRESH, 
                        FREQUENCY_THRESHOLD=FREQUENCY_THRESHOLD, GENETIC_CODE=GENETIC_CODE,
                        COUNT_THRESHOLD=COUNT_THRESHOLD, LIMIT_TO_SUSPECTED_SPECIES=LIMIT_TO_SUSPECTED_SPECIES)

            q.put((args, setting, runid ))

        else:
            res['error'] = True
            BadRunList.put(runid)
    else :
        print 'Error is from request'
        res['error'] = True

    return json.dumps(res)


def save_uploaded_file(folder, filename, file, chunk_size=1000):
    """ save file in a directory"""
    # Using chung reduce the data size
    # saving memory
    with open(os.path.join(folder,filename), 'wb') as FILE:
        while True:
            data = file.read(chunk_size)
            if not data : break
            FILE.write(data)

def index(environ, start_response):
    """Return the index file mounted at "/" and display"""
    start_response('200 OK', [('Content-Type', 'text/html')])
    return [open(os.path.join(WEB_APP_BASE_PATH, "webplugin/index.html"), 'r').read()]

def not_found(environ, start_response):
    """Called if no URL matches."""
    start_response('404 NOT FOUND', [('Content-Type', 'text/plain')])
    return ['Not Found']


def uploadfile(environ, start_response):
    """ Upload file here"""
    start_response('200 OK', [('Content-Type', 'application/json')])
    result = {}
    t = time.time() # timestamp to find url relative to the current run
    try : 
        post_env = environ.copy()
        post = cgi.FieldStorage(
            fp=post_env['wsgi.input'],
            environ=post_env,
            keep_blank_values=True
        )
        folder = os.path.join(WEB_APP_BASE_PATH, 'input', str(t))
        if not os.path.exists(folder):
            os.makedirs(folder)
        for f in post.keys():
            c_file = post[f].file
            if c_file:
                save_uploaded_file(folder, f, c_file)
                # returning all the path is a huge security risk
                # let's just return the relative path
                result[f] = os.path.join(str(t), f)

    except Exception as e:
        print e
        result['error'] = "Problem with uploading your file. Please try again.\nContact the admin if the problem persist"

    result['time'] = str(t) 
    return json.dumps(result)


def track(environ, start_response):
    # expected arguments from the URL (POST or GET method)

    queries = {}
    
    if environ['REQUEST_METHOD'].upper() == 'GET' and  environ['QUERY_STRING']:
        queries = cgi.parse_qs(environ['QUERY_STRING'])
    elif environ['REQUEST_METHOD'].upper() == 'POST' and environ['wsgi.input']:
        queries = cgi.parse_qs(environ['wsgi.input'].read())

    param_seqid = queries.get("seqtype", [None])[0]
    treeid = queries.get("treeid", [None])[0]

    param_browser_id_start = queries.get("browser_id_start", [None])[0]
    if None in set([param_seqid]):
        return ["Error, Not enough params"]

    if not treeid:
        treeid = md5(str(time.time())).hexdigest()

    else:
        sptree = PhyloTree(TMP+TREE, alignment=TMP+ALIGNMENT[param_seqid])
        
        selected_align = AlignIO.read(TMP+ALIGNMENT[param_seqid], 'fasta', alphabet=Alphabet.Gapped(IUPAC.protein))
        summary_info = AlignInfo.SummaryInfo(selected_align)        
        total_ic_content = summary_info.information_content()
        global IC_CONTENT
        IC_CONTENT = summary_info.ic_vector.values()
        threshold = float(max(IC_CONTENT)* settings.IC_INFO_THRESHOLD)
        #print threshold, max(IC_CONTENT)

        global ts
        ts = TreeStyle()
        ts.branch_vertical_margin = 5
        ts.show_scale = False
        ic_plot = faces.SequencePlotFace(IC_CONTENT, hlines=[threshold], hlines_col=['red'], ylim=(int(min(IC_CONTENT)-0.5), int(max(IC_CONTENT)+0.5)), fsize=10, col_width=14, header="Information Content", kind='bar', ylabel="ic")
        ts.aligned_header.add_face(ic_plot, 1)            
        application.set_tree_style(ts)

        newick = sptree.write(features=[])
        result =""" <script>$("#treenewick").css("visibility", "visible");
                        $("#treenewick").val("%s");
                        draw_tree(%s, "%s", "#img1");
                        $("#treenewick").attr('disabled', true);

                        </script>"""%(newick, treeid, newick)

        return [result]


def checkrun(environ, start_response):
    """This serve to check each run from the link"""
    queries = {}
    if environ['REQUEST_METHOD'].upper() == 'GET' and  environ['QUERY_STRING']:
        queries = cgi.parse_qs(environ['QUERY_STRING'])
        print queries
    elif environ['REQUEST_METHOD'].upper() == 'POST' and environ['wsgi.input']:
        queries = cgi.parse_qs(environ['wsgi.input'].read())


    runid = queries.get('id', [None])[0]
    global ThreadParam
    valid_run = ThreadParam.RES
    exec_run = ThreadParam.EXECUTING
    failed_run = ThreadParam.ERRORS

    def runnotfound(runid):
        print runid
        print q.queue
        print exec_run
        print valid_run
        print failed_run
        return not (runid and (runid in q.queue or runid in valid_run.keys() or runid in failed_run or runid in exec_run))

    if runnotfound(runid):
        return not_found(environ, start_response)

    elif runid in valid_run.keys():
        start_response('202 OK', [('Content-Type', 'text/plain')])
        return valid_run[runid].get_result()
    
    elif runid in exec_run:
        start_response('202 OK', [('Content-Type', 'text/plain')])
        return ['Refresh this page later']
    else:
        start_response('202 OK', [('Content-Type', 'text/plain')])
        return ['There was an error with your request']

    
# map urls to functions
urls = [
    (r'^$', index),
    (r'/track(.*)$', track),
    (r'/uploadfile(.*)$', uploadfile),
    (r'/calltracker(.*)$', coretracker_run),
    (r'/results/run(.*)$', checkrun),
]

def application(environ, start_response):
    asked_method = environ.get('PATH_INFO', '').lstrip('/')
    print "REALLY"
    for regex, callback in urls:
        match = re.search(regex, asked_method)
        if match is not None:
            return callback(environ, start_response)

    # maybe use json output here
    return not_found(environ, start_response)

        

# ==============================================================================
# TREE LOADER
#
# This is my tree loading functions. I want the WebTreeApplication to
# use this method to load the trees
# ==============================================================================


def phyloloader(tree):
    """ This is function is used to load trees within the
    WebTreeApplication object. """
    t = PhyloTree(tree)
    return t

# ==============================================================================
# CUSTOM LAYOUTS
#
# This are my layout functions. I want the WebTreeApplication to use
# them for rendering trees
# ==============================================================================

LEAVE_FACES = [] # Global var that stores the faces that are rendered
                 # by the layout function

def main_layout(node):
    ''' Main layout function. It controls what is shown in tree
    images. '''

    # Add faces to leaf nodes. This allows me to add the faces from
    # the global variable LEAVE_FACES, which is set by the application
    # controler according to the arguments passed through the URL.
    global IC_CONTENT

    if not (IC_CONTENT):
        print len(IC_CONTENT)

        selected_align = AlignIO.read(TMP+ALIGNMENT[param_seqid], 'fasta', alphabet=Alphabet.Gapped(IUPAC.protein))
        summary_info = AlignInfo.SummaryInfo(selected_align)        
        total_ic_content = summary_info.information_content()
        IC_CONTENT = summary_info.ic_vector.values()

    global ts
    ic_plot = faces.SequencePlotFace(IC_CONTENT, fsize=10, col_width=14, header="Information Content", kind='bar', ylabel="ic")
    ts.aligned_header.add_face(ic_plot, 1)
 
    
    if node.is_leaf():
        for f, fkey, pos in LEAVE_FACES:
            if hasattr(node, fkey) and fkey != 'name' :
                # we check if the face is already added before 
                print fkey
                if(fkey=='name'):
                    added = getattr(node, '_temp_faces', None)
                    print getattr(added, 'branch-right')
                faces.add_face_to_node(f, node, column=pos, position="branch-right")
    else:
        # Add special faces on collapsed nodes
        if hasattr(node, "hide") and int(node.hide) == 1:
            print "trouble"

            node.img_style["draw_descendants"]= False
            collapsed_face = faces.TextFace(\
                " %s collapsed leaves." %len(node), \
                    fsize=10, fgcolor="#444", ftype="Arial")
            faces.add_face_to_node(collapsed_face, node, 0)
        else:
            node.img_style["draw_descendants"] = True

    
    # Set node aspect. This controls which node features are used to
    # control the style of the tree. You can add or modify this
    # features, as well as their behaviour
    if node.is_leaf():
        node.img_style["shape"] = "square"
        node.img_style["size"] = 4
    elif node.is_root():
        node.img_style["shape"] = "circle"
        node.img_style["size"] = 8
    else:
        node.img_style["size"] = 8
        node.img_style["shape"] = "sphere"

    
    # Parse node features features and conver them into styles. This
    # must be done like this, since current ete version does not allow
    # modifying style outside the layout function.
    """
    if hasattr(node, "bsize"):
        node.img_style["size"]= int(node.bsize)

    if hasattr(node, "shape"):
        node.img_style["shape"]= node.shape

    if hasattr(node, "bgcolor"):
        node.img_style["bgcolor"]= node.bgcolor

    if hasattr(node, "fgcolor"):
        node.img_style["fgcolor"]= node.fgcolor
    """

# ==============================================================================
# Checker function definitions:
#
# All checker actions must receive a node instance as unique argument
# and return True (node passes the filters) or False (node does not
# passes the filters).
#
# ==============================================================================

can_expand = lambda node: not node.is_leaf() and (hasattr(node, "hide") and node.hide==True)
can_collapse = lambda node: not node.is_leaf() and (not hasattr(node, "hide") or node.hide==False)
is_leaf = lambda node: node.is_leaf()
is_not_leaf = lambda node: not node.is_leaf()

# ==============================================================================
# Handler function definitions:
#
# All action handler functions must receive a node instance as unique
# argument. Returns are ignored.
#
# Note that there is a special action handler designed for searches
# within the tree. Handler receives node and searched term.
#
# ==============================================================================

def collapse(node):
    node.add_feature("hide", 1)
    node.add_feature("bsize", 25)
    node.add_feature("shape", "sphere")
    node.add_feature("fgcolor", "#bbbbbb")

def expand(node):
    try:
        node.del_feature("hide")
        node.del_feature("bsize")
        node.del_feature("shape")
        node.del_feature("fgcolor")
    except (KeyError, AttributeError):
        pass

def swap_branches(node):
    node.children.reverse()

def set_red(node):
    node.add_feature("fgcolor", "#ff0000")
    node.add_feature("bsize", 40)
    node.add_feature("shape", "sphere")

def set_bg(node):
    node.add_feature("bgcolor", "#CEDBC4")

def set_as_root(node):
    node.get_tree_root().set_outgroup(node)

def phylomedb_clean_layout(node):
    phylomedb_layout(node)
    node.img_style["size"]=0

def search_by_feature(tree, search_term):
    ''' Special action '''
    attr, term = search_term.split("::")
    if not term:
        return None
    elif attr == "clean" and term == "clean":
        for n in tree.traverse():
            try:
                n.del_feature("bsize")
                n.del_feature("shape")
                n.del_feature("fgcolor")
            except:
                pass
    else:
        for n in tree.traverse():
            if hasattr(n, attr) and \
                    re.search(term,  str(getattr(n, attr)), re.IGNORECASE):
                n.add_feature("bsize", 16)
                n.add_feature("shape", "sphere")
                n.add_feature("fgcolor", "#BB8C2B")


# ==============================================================================
# HTML generators
#
# Actions will be automatically added to the popup menus and attached
# to the action handler function. However, if just want to add
# informative items to the popup menu or external actions not
# associated to any handler, you can overwrite the default html
# generator of each action.
#
# html generators receive all information attached to the node and
# action in 5 arguments:
#
# * aindex: index of the action associated to this html generator
#
# * nodeid: id of the node to which action is attached
#
# * treeid: id of the tree in which node is present
#
# * text: the text string associated to the element that raised the
# action (only applicable to text faces actions)
#
# * node: node instance in which action will be executed.
#
#
# Html generator should return a text string encoding a html list
# item:
#
# Example: return "<li> my text </li>"
#
# ==============================================================================


def branch_info(aindex, nodeid, treeid, text, node):
    ''' It shows some info of the node in the popup menu '''
    return """
           <li style="background:#eee; font-size:8pt;">
           <div style="text-align:left;font-weight:bold;">
            NODE ACTIONS
           </div>
            (<b>Branch: </b>%0.3f <b>Support:</b> %0.3f)<br>
           </li>"""  %\
        (node.dist, node.support)

def search_in_ensmbl(aindex, nodeid, treeid, text, node):
    return '''<li>
              <a target="_blank" href="http://www.ensembl.org/common/Search/Results?species=all;idx=;q=%s">
              <img src=""> Search in ensembl: %s >
              </a>
              </li> ''' %\
            (node.name, node.name)

def external_links_divider(aindex, nodeid, treeid, text, node):
    ''' Used to show a separator in the popup menu'''
    if node.is_leaf():
        return """<li
        style="background:#eee;font-size:8pt;"><b>External
        links</b></li>"""
    else:
        return ""

def topology_action_divider(aindex, nodeid, treeid, text, node):
    return """<li style="background:#eee;"><b>Tree node actions</b></li>"""

# ==============================================================================
# TREE RENDERER
#
# By default, ETE will render the tree as a png image and will return
# a simplistic HTML code to show the image and interact with
# it. However, it is possible to wrap such functionality to preprocess
# trees in a particular way, read extra parameters from the URL query
# and/or produce enriched HTML applications.
#
# Tree renderer wrappers receive the tree object, its id, and the WSGI
# application object. They MUST call the
# application._get_tree_img(tree) method and return the a HTML
# string.
#
# A simplistic wrapper that emulates the default WebTreeApplication
# behaviour would be:
#
# def tree_renderer(tree, treeid, application):
#    html = application._get_tree_img(treeid = treeid)
#    return html
#
# ==============================================================================
def tree_renderer(tree, treeid, application):
    # The following part controls the features that are attched to
    # leaf nodes and that will be shown in the tree image. Node styles
    # are set it here, and faces are also created. The idea is the
    # following: user can pass feature names using the URL argument
    # "tree_features". If the feature is handled by our function and
    # it is available in nodes, a face will be created and added to
    # the global variable LEAVE_FACES. Remember that our layout
    # function uses such variable to add faces to nodes during
    # rendering.

    # Extracts from URL query the features that must be drawn in the tree 
    asked_features = application.queries.get("show_features", ["name"])[0].split(",")

    def update_features_avail(feature_key, name, col, fsize, fcolor, prefix, suffix):
        text_features_avail.setdefault(feature_key, [name, 0, col, fsize, fcolor, prefix, suffix])
        text_features_avail[feature_key][1] += 1

    tree.add_feature("fgcolor", "#000000")
    tree.add_feature("shape", "sphere")
    tree.add_feature("bsize", "8")
    tree.dist = 0

    # This are the features that I wanto to convert into image
    # faces. I use an automatic function to do this. Each element in
    # the dictionary is a list that contains the information about how
    # to create a textFace with the feature.
    leaves = tree.get_leaves()
    formated_features =  {
        # feature_name: ["Description", face column position, text_size, color, text_prefix, text_suffix ]
        "name": ["Leaf name", len(leaves), 0, 12, "#000000", "", ""],
        "spname": ["Species name", len(leaves), 1, 12, "#f00000", " Species:", ""],
        }

    # populates the global LEAVE_FACES variable
    global LEAVE_FACES
    LEAVE_FACES = []
    unknown_faces_pos = 2
    for fkey in asked_features:
        if fkey in formated_features:
            name, total, pos, size, color, prefix, suffix = formated_features[fkey]
            f = faces.AttrFace(fkey, ftype="Arial", fsize=size, fgcolor=color, text_prefix=prefix, text_suffix=suffix)
            LEAVE_FACES.append([f, fkey, pos])
        else:
            # If the feature has no associated format, let's add something standard
            prefix = " %s: " %fkey
            suffix = ","
            if fkey == asked_features[-1]:
                suffix=""
            f = faces.AttrFace(fkey, ftype="Arial", fsize=10, fgcolor="#666666", text_prefix=prefix, text_suffix=suffix)
            LEAVE_FACES.append([f, fkey, unknown_faces_pos])
            unknown_faces_pos += 1

    text_features_avail = {}
    for l in leaves:
        for f in l.features: 
            if not f.startswith("_"):
                text_features_avail.setdefault(f, 0)
                text_features_avail[f] = text_features_avail[f] + 1 

    html_features = '''
      <div id="tree_features_box">
      <div class="tree_box_header">Available tree features
      <img src=""/CoreTracker/webplugin/images/close.png" onclick='$(this).closest("#tree_features_box").hide();'>
      </div>
      <form id="form_tree_features" action='javascript: set_tree_features("", "", "");'>
      '''

    for fkey, counter in text_features_avail.iteritems():
        if fkey in asked_features:
            tag = "CHECKED"
        else:
            tag = ""

        fname = formated_features.get(fkey, [fkey])[0]

        html_features += '<INPUT NAME="tree_feature_selector" TYPE=CHECKBOX %s VALUE="%s">%s (%s/%s leaves)</input><br> ' %\
            (tag, fkey, fname, counter, len(leaves))


    html_features += """<input type="submit" value="Refresh" 
                        onclick='javascript:
                                // This piece of js code extracts the checked features from menu and redraw the tree sending such information
                                var allVals = [];
                                $(this).parent().children("input[name=tree_feature_selector]").each(function(){
                                if ($(this).is(":checked")){
                                    allVals.push($(this).val());
                                }});
                                draw_tree("%s", "", "#img1", {"show_features": allVals.join(",")} );'
                       > 
                       </form></div>""" %(treeid)

    features_button = """
     <li><a href="#" onclick='show_box(event, $(this).closest("#tree_panel").children("#tree_features_box"));'>
     <img width=16 height=16 src="/CoreTracker/webplugin/images/icon_tools.png" alt="Select Tree features">
     </a></li>"""

    download_button = """
     <li><a href="/CoreTracker/tmp/%s.png" target="_blank">
     <img width=16 height=16 src="/CoreTracker/webplugin/images/icon_attachment.png" alt="Download tree image">
     </a></li>""" %(treeid)

    search_button = """
      <li><a href="#" onclick='javascript:
          var box = $(this).closest("#tree_panel").children("#search_in_tree_box");
          show_box(event, box); '>
      <img width=16 height=16 src="/CoreTracker/webplugin/images/icon_search.png" alt="Search in tree">
      </a></li>"""

    clean_search_button = """
      <li><a href="#" onclick='run_action("%s", "", %s, "clean::clean");'>
      <img width=16 height=16 src="/CoreTracker/webplugin/images/icon_cancel_search.png" alt="Clear search results">
      </a></li>""" %\
        (treeid, 0)

    buttons = '<div id="ete_tree_buttons">' +\
        features_button + search_button + clean_search_button + download_button +\
        '</div>'

    search_select = '<select id="ete_search_target">'

    # Name first
    if 'name' in text_features_avail:
        search_select += '<option value="name">name</option>'
        del text_features_avail['name']

    for fkey in text_features_avail:
        search_select += '<option value="%s">%s</option>' %(fkey,fkey)
    search_select += '</select>'


    main_search = '''
    <form onsubmit='javascript:
                     search_in_tree("%s", "%s",
                                    $(this).closest("form").children("#ete_search_term").val(),
                                    $(this).closest("form").children("#ete_search_target").val());'
          action="javascript:void(0);">
          <input id="ete_search_term" type="text" value="Search"
          style="font-style:italic;color:grey;"
          onfocus="if(this.value == 'Search') { this.value = ''; this.style.color='black'; this.style.fontStyle='normal'}"'>%s</input>
     </form>
     '''%\
            (treeid, 0, search_select)


    search_form = """
     <div id="search_in_tree_box">
     <div class="tree_box_header"> Search in Tree
     <img src="/CoreTracker/webplugin/images/close.png" onclick='$(this).closest("#search_in_tree_box").hide();'>
     </div>
     <form onsubmit='javascript:
                     search_in_tree("%s", "%s",
                                    $(this).closest("form").children("#ete_search_term").val(),
                                    $(this).closest("form").children("#ete_search_target").val());'
          action="javascript:void(0);">
     %s
     <input id="ete_search_term" type="text" value=""'></input>
     <br><i>(Searches are not case sensitive and accept Perl regular expressions)</i>
     <br>
     </form>
     <i> (Press ENTER to initiate the search)</i>
     </div>
     """ %\
        (treeid, 0, search_select) # 0 is the action index associated
                                   # to the search functionality. This
                                   # means that this action is the
                                   # first to be registered in WebApplication.
    
    buttons = '<div id="ete_tree_buttons">' +\
            main_search + features_button + clean_search_button + download_button +\
            '</div>'

    tree_panel_html = '<div id="tree_panel">' + search_form + html_features + buttons + '</div>'

    # Now we render the tree into image and get the HTML that handles it
    tree_html = application._get_tree_img(treeid = treeid)

    # Let's return enriched HTML
    return tree_panel_html + tree_html

# ==============================================================================
#
# Main WSGI Application
#
# ==============================================================================

# Create a basic ETE WebTreeApplication
application = WebTreeApplication(application)

# this is another thread to deleted all file older than 3days
del_thread = RunDelThread(del_run, EXPIRE_TIME, CHECK_TIME, BadRunList)
del_thread.setDaemon(True)
del_thread.start()
del_thread.name = 'Request_run_gc'

# we start 5 thread to run request
for i in xrange(THREAD_MAX):
    print "this is hella strange ", i
    t = CoreTrackerThread(call_coretracker, q)
    t.setDaemon(True)
    t.start()
    t.name = "Coretracker_exec_%d"%i

# Set your temporal dir to allow web user to generate files. This two
# paths should point to the same place, one using the absolute path in
# your system, and the other the URL to access the same
# directory. Note that the referred directory must be writable by the
# webserver.
#application.CONFIG["temp_dir"] = "/home/services/web/etetoolkit.org/webplugin/tmp/"
application.CONFIG["temp_dir"] = "/home/manu/html/application/tmp"
application.CONFIG["temp_url"] = "./../tmp" # Relative to web site Document Root

# Set the DISPLAY port that ETE should use to draw pictures. You will
# need a X server installed in your server and allow webserver user to
# access the display. If the X server is started by a different user
# and www-data (usally the apache user) cannot access display, try
# modifiying DISPLAY permisions by executing "xhost +"
application.CONFIG["DISPLAY"] = ":0" # This is the most common
                                     # configuration

 # start virtual display
#xvfb=Display(visible=0, size=(1024, 768)).start()
#application.CONFIG["DISPLAY"] = str(xvfb.new_display_var)

# Lets now apply our custom tree loader function to the main
# application 
application.set_tree_loader(phyloloader)

# And our layout as the default one to render trees

#ts.layout_fn.append(main_layout)
#application.set_tree_style(ts)

#application.set_default_layout_fn(main_layout)
# application.set_tree_size(None, None)
# I want to make up how tree image in shown using a custrom tree
# renderer that adds much more HTML code
application.set_external_tree_renderer(tree_renderer)


# ==============================================================================
# ADD CUSTOM ACTIONS TO THE APPLICATION
#
# The function "register_action" allows to attach functionality to
# nodes in the image. All registered accions will be shown in the
# popup menu bound to the nodes and faces in the web image.
#
#
# register_action(action_name, target_type=["node"|"face"|"layout"|"search"], \
#                 action_handler, action_checker, html_generator_handler)
#
# When the Application is executed it will read your registered
# acctions and will do the following:
#
# 1. Load the tree and get the image map
#
# 2. For each node and face in the tree, it will browse all registered
# actions and will run the action_checker function to determine if the
# action must be activated for such node or face
#
# 3. If action_checker(node) returns True, the action will be attached
# to the context menu of that specific node or face, otherwise it will
# be hidden.
#
# 4. When a click is done on a specific node, popup menus will be
# built using their active actions. For this, ETE will use the
# html_generator function associated to each function if
# available. Otherwise, a popup entry will be added automatically.
#
# 5. When a certain action is pressed in the popup menus, the
# action_handler function attached to the action will be executed over
# its corresponding node, and the tree image will be refreshed.
#
# Special values:
#
#  action_checker = None : It will be interpreted as "Show allways"
#  html_generator = None : Autogenerate html and link to action
#  action_handler = None : Action will be ignored
#
# ==============================================================================

# We first register the special action "search" which is attached to
# our custom search function.
application.register_action("", "search", search_by_feature, None, None)

# Node manipulation options (bound to node items and all their faces)
application.register_action("branch_info", "node", None, None, branch_info)
application.register_action("<b>Collapse</b>", "node", collapse, can_collapse, None)
application.register_action("Expand", "node", expand, can_expand, None)
application.register_action("Highlight background", "node", set_bg, None, None)
application.register_action("Set as root", "node", set_as_root, None, None)
application.register_action("Swap children", "node", swap_branches, is_not_leaf, None)

# Actions attached to node's content (shown as text faces)
#application.register_action("divider", "face", None, None, external_links_divider)
#xvfb.stop()