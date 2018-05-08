from template import BasicTemplate
from weasyprint import HTML
import os

curdir = os.path.dirname(os.path.realpath(__file__))
template_dir = os.path.join(curdir, 'templates')
templates = [os.path.join(template_dir, "default.html"),
             os.path.join(template_dir, "rea_improve.html")]
default_style = os.path.join(template_dir, "style.css")


def export_from_html(template_vars, output, csslist=[default_style], base_url=None, template=0):
    """export from html to pdf using weasyprint"""
    if not isinstance(csslist, list):
        raise ValueError("Expect list for csslist, got " + type(csslist) + "!")
    template = templates[template]
    tmp = BasicTemplate(template)
    inputhtml = tmp.render(template_vars)
    HTML(string=inputhtml.encode('utf-8'), base_url=base_url).write_pdf(output, stylesheets=csslist)
    return output
