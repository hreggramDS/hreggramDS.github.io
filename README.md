# Greg Hamilton's Data Science Portfolio

[![GitHub Pages](https://img.shields.io/badge/GitHub-Pages-blue)](https://hreggramds.github.io)

A Jekyll-based portfolio site showcasing data science projects, research, and technical writing.

## About This Site

This site serves as a professional portfolio and research hub, featuring:
- **Projects**: Detailed documentation of data science and research projects
- **Pipelines**: Technical documentation for reproducible workflows
- **Blog**: Technical tutorials and insights
- **Publications**: Academic papers and presentations

## Local Development

### Prerequisites
- Ruby 2.7 or higher
- Bundler

### Setup

1. Clone the repository:
```bash
git clone https://github.com/hreggramDS/hreggramDS.github.io.git
cd hreggramDS.github.io
```

2. Install dependencies:
```bash
bundle install
```

3. Run the local server:
```bash
bundle exec jekyll serve
```

4. View the site at `http://localhost:4000`

### Building the Site

```bash
bundle exec jekyll build
```

The built site will be in the `_site/` directory.

## Site Structure

```
hreggramDS.github.io/
├── _config.yml           # Site configuration
├── _data/                # Structured data (YAML)
├── _includes/            # Reusable components
├── _layouts/             # Page templates
├── _posts/               # Blog posts (YYYY-MM-DD-title.md)
├── _projects/            # Project collection
├── _pipelines/           # Pipeline documentation
├── assets/               # Static assets (CSS, JS, images)
├── pages/                # Static pages
├── index.md              # Landing page
└── SITE_DESIGN.md        # Site design documentation
```

## Adding Content

### New Blog Post

Create a file in `_posts/` with the format `YYYY-MM-DD-title.md`:

```markdown
---
layout: post
title: "Your Post Title"
date: YYYY-MM-DD
categories: ["category1", "category2"]
tags: ["tag1", "tag2"]
author: "Greg Hamilton"
excerpt: "Brief description"
---

Your content here...
```

### New Project

Create a file in `_projects/` with your project name:

```markdown
---
layout: project
title: "Project Title"
date: YYYY-MM-DD
status: "Completed"
domain: ["Machine Learning"]
technologies: ["Python", "scikit-learn"]
github: "https://github.com/username/repo"
description: "Brief description"
featured: true
---

Your project details...
```

### New Pipeline

Create a file in `_pipelines/`:

```markdown
---
layout: pipeline
title: "Pipeline Name"
date: YYYY-MM-DD
version: "1.0"
domain: "Data Engineering"
technologies: ["Python", "Airflow"]
repository: "https://github.com/username/repo"
---

Your pipeline documentation...
```

## Deployment

This site is automatically deployed to GitHub Pages when changes are pushed to the `main` branch.

## Technologies

- **Jekyll**: Static site generator
- **GitHub Pages**: Hosting
- **Markdown**: Content format
- **Liquid**: Templating language
- **Jekyll Theme Slate**: Base theme with custom enhancements

## Documentation

See [SITE_DESIGN.md](SITE_DESIGN.md) for comprehensive design documentation and implementation details.

## License

This project is licensed under the terms specified in [LICENSE](LICENSE).

## Contact

- **GitHub**: [@hreggramDS](https://github.com/hreggramDS)
- **Website**: [hreggramds.github.io](https://hreggramds.github.io)

---

Built with ❤️ using Jekyll and GitHub Pages
