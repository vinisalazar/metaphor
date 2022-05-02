"""
Generate docs for the config module.
"""
default_config = "metaphor/config/default_value-config.yaml"
docs_page = "docs/source/main/settings.md"

final_md = """
# Configuration

This page explains each value of Metaphor's config settings, that is, the values defined in the config YAML file.
\n
"""


def parse_module_settings(string):
    module_settings = [line.strip() for line in string.split("\n\n") if line]
    module_settings = [line for line in module_settings if not line.startswith("#")]
    return module_settings


def parse_config_property(string):
    properties = [line.strip() for line in string.split("\n")]
    fmt_properties = []
    for prop in properties:
        try:
            prop_name, remainder = [i.strip() for i in prop.split(":", 1)]
        except:
            breakpoint()
        try:
            if "#" in remainder:
                prop_default_value, prop_comment = [
                    i.strip() for i in remainder.split("#")
                ]
            else:
                prop_default_value, prop_comment = remainder.strip(), ""
        except:
            breakpoint()
        fmt_properties.append((prop_name, prop_default_value, prop_comment))

    return fmt_properties


with open(default_config) as f:
    text = f.read()
    text = text.split("###############################################################")
    text = [tuple(text[line : line + 2]) for line in range(1, len(text), 2)]


for line in text:
    module, settings = line
    module = [i.replace("# ", "") for i in module.split("\n") if i]
    module_name, *lines = module

    final_md += f"**{module_name}**\n\n"
    for line in lines:
        final_md += f"{line}\n\n"

    module_settings = parse_module_settings(settings)
    config_properties = [
        parse_config_property(config_property) for config_property in module_settings
    ]
    for item in config_properties:
        tab = False
        for name, default_value, comment in item:
            name_fmt = f"**`{name}:`**"
            if tab:
                name_fmt = f"&nbsp;&nbsp;&nbsp;{name_fmt}"
            final_md += name_fmt
            if default_value:
                final_md += f" `{default_value}`  "
            else:
                final_md += "  "
                tab = True
            if comment:
                final_md += f" {comment}  "
            else:
                final_md += "  "
            final_md += "\n"

        final_md += "\n\n"


with open(docs_page, "w") as f:
    f.write(final_md)
