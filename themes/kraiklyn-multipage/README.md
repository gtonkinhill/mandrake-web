Kraiklyn-multipage
=====

A simple [Hugo](https://gohugo.io/) theme for documentation. Based on the Kraiklyn one page theme

# Features
 - unlimited menu level
 - mobile friendly
 - customize website logo
 - add links to the sidebar

# Installation
Clone the repository to your site’s themes directory.

# Usage
All content is rendered on it's own page. Content is ordered by Weight.


## Shortcodes

### block
Create notes, tips and other blocks on the page
```markdown
{{% block note %}}
By default only ports 22, 80 and 443 are open
{{% /block %}}
```

Available types: `note`, `tip`, `warn`, `info`


## Customizing sidebar

### Changing logo
Replace logo by creating `layouts/partials/logo.html` file

### Adding menu entries to the external links section
Customize the name of the section by adding to the `config.toml`
```toml
[params]
externalTitle = "Surfly docs"
```

Add new entries:
```toml
[[menu.shortcuts]]
name = "Javascript API"
url = "https://docs.surfly.com/javascript.html"
weight = 20
```

### Changing color
Customize the color of the sidebar by adding to the `config.toml`
```toml
[params]
sidebarColor = "green"
```

Available values : default, green, purple, pink, red, cyan, blue, grey, orange.

## Add favicon
Put `favicon.ico` inside `static` folder


## Add custom CSS
You can add your custom CSS files with the `customCss` parameter of the configuration file.

```toml
customCss = ["css/custom.css", "css/custom2.css"]
```

Just put your files in `static/css` directory.
