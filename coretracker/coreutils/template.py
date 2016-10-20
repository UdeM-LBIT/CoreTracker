from output import Output


class BasicTemplate(object):
    """Really basic templating system
    This is just to prevent too much dependency when installing
    """

    def __init__(self, template):
        self.template = ""
        self._get_template(template)

    def _get_template(self, template):
        """Get template from file"""
        with open(template, 'r') as INPUT:
            # read is unreliable
            self.template = "\n".join(INPUT.readlines())

    def render(self, varlist):
        """Render and return the string"""
        return self.template.format(**varlist)

    def render_save(self, varlist, outputfile=None):
        """Render and save to file"""
        rendered_str = self.render(varlist)
        out = Output(outputfile)
        out.write(rendered_str)
        out.close()
        return outputfile
