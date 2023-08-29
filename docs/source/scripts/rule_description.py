"""
This script auto-generates the rule description for each rule in the workflow.

This was taken from the seq2science repository, under an MIT licence.

Link to repository: https://github.com/vanheeringen-lab/seq2science
"""
short_license = """
---

**Disclaimer**

This page was generated with a script adapted from the 
[seq2science repository](https://github.com/vanheeringen-lab/seq2science).

MIT License

Copyright (c) 2019 Maarten-vd-Sande (vanheeringen-lab)

For the full license, please see the
[script source code](https://github.com/vinisalazar/metaphor/blob/main/docs/scripts/rule_description.py).
"""

license = """
This page was generated with a script adapted from the 
[seq2science repository](https://github.com/vanheeringen-lab/seq2science). Please find the license attached.\n

**Licence for rule_description.py script:**

MIT License

Copyright (c) 2019 Maarten-vd-Sande (vanheeringen-lab)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

"""

import re

final_md = """
# Reference

This is an automatically generated list of all supported rules, their docstrings, and command. At the start of each 
workflow run a list is printed of which rules will be run. And while the workflow is running it prints which rules are
being started and finished. This page is here to give an explanation to the user about what each rule does, and for
developers to find what is, and isn't yet supported. Not all Metaphor rules are listed here, only the ones with a
`shell` directive. Rules with `script` or `wrapper` directives are not included. To see all rules in Metaphor, 
please refer to the [workflow source code](https://github.com/vinisalazar/metaphor/tree/main/workflow).
"""

path = "metaphor/workflow/rules/"
docs_page = "docs/source/main/reference.md"


def get_dirty_docstrings(string):
    splitter = re.compile(
        r"(rule|checkpoint) (.*):[\s\S]*?\"\"\"([\s\S]*?)\"\"\"", re.MULTILINE
    )
    docstrings = {}
    for match in splitter.finditer(string):
        docstrings[match.group(2)] = match.group(3)
    return docstrings


def cleanup_docstring(dirty):
    clean = {}
    for rule, docstring in dirty.items():
        if len(firstline := docstring.split("\n")) > 1:
            firstline = firstline[1]
        else:
            firstline = firstline[0]
        indentation = len(firstline) - len(firstline.lstrip())
        docstring = docstring.replace(" " * indentation, "")
        docstring = docstring.replace(" " * (indentation - 4), "")
        docstring = docstring.strip("\n")
        docstring = docstring.replace("#", "\#")
        clean[rule] = docstring

    return clean


def get_dirty_shell(string):
    splitter = re.compile(
        r"(rule|checkpoint) (.*):[\s\S]*?(shell|script|wrapper):[\s\S]*?\"\"\"\\?([\s\S]*?)\"\"\"",
        re.MULTILINE,
    )
    shell_cmds = {}
    for substring in string.split("\n\n\n"):
        for match in splitter.finditer(substring):
            shell_cmds[match.group(2)] = match.group(4)
    return shell_cmds


def cleanup_shell(dirty):
    clean = {}
    for rule, shell in dirty.items():
        if len(firstline := shell.split("\n")) > 1:
            firstline = firstline[1]
        else:
            firstline = firstline[0]
        indentation = len(firstline) - len(firstline.lstrip())
        shellstring = "\n".join(
            [
                shell_line.replace(" " * indentation, "", 1)
                for shell_line in shell.split("\n")
            ]
        )
        shellstring = shellstring.strip("\n")
        # shellstring = shellstring.replace("#", "\#")
        clean[rule] = shellstring

    return clean


all_rules_doc = {}
all_rules_shell = {}

rule_orders = {
    "qc": None,
    "assembly": None,
    "annotation": None,
    "mapping": None,
    "binning": None,
    # "postprocessing": None,
}

for rules_file in rule_orders.keys():
    with open(path + rules_file + ".smk", "r") as file:
        text = file.read()
    shell_cmd = cleanup_shell(get_dirty_shell(text))
    all_rules_shell.update(shell_cmd)

    docstrings = cleanup_docstring(get_dirty_docstrings(text))
    all_rules_doc.update(docstrings)
    rule_orders[rules_file] = list(docstrings.keys())

    final_md += f"## {rules_file}.smk\n"

    for rule in rule_orders[rules_file]:
        print(f"Writing docs for rule '{rule}'.")
        docstring = all_rules_doc[rule]
        final_md += f"**{rule}**\n\n"
        if "{" not in docstring:
            final_md += f"{docstring}\n\n"
        if rule in all_rules_shell:
            final_md += "```\n"
            final_md += f"{all_rules_shell[rule]}\n"
            final_md += "```\n\n"

final_md += short_license

with open(docs_page, "w") as text_file:
    text_file.write(final_md)
