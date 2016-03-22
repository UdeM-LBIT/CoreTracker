from PyPDF2 import PdfFileReader, PdfFileWriter
from template import BasicTemplate
from weasyprint import HTML
import os

template_dir =  os.path.abspath('../templates')
templates = [os.path.join(template_dir, "default.html"),
                os.path.join(template_dir, "rea_improve.html")]
default_style = os.path.join(template_dir, "style.css")

def export_from_html(template_vars, output, csslist=[default_style], base_url=None, template=0):
    """export from html to pdf using weasyprint"""
    if not isinstance(csslist, list):
        raise ValueError("Expect list for csslist, got "+type(csslist)+"!")
    template = templates[template]
    tmp = BasicTemplate(template)
    inputhtml = tmp.render(template_vars)
    HTML(string=inputhtml, base_url=base_url).write_pdf(output, stylesheets=csslist)
    return output

def concat_pdf(outfile, *args):
    """Concat multiple pdf into one"""
    output = PdfFileWriter()
    npdf = len(args)
    if npdf < 1:
        return
    else :
        args = [PdfFileReader(arg) for arg in args]
        arg1 = args[0]
        dim = arg1.getPage(0).mediaBox
        for arg in args:
            for p in arg.pages:
                p.scaleTo(dim[2], dim[3])
                output.addPage(p)

    with open(outfile, "wb") as OUTSTREAM:
        output.write(OUTSTREAM)
