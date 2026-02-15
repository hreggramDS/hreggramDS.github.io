---
layout: page
title: Projects
description: Portfolio of data science projects, research, and technical work
permalink: /pages/projects/
---

# Projects Portfolio

Explore my data science and research projects across various domains including bioinformatics, machine learning, and scientific computing.

## Featured Projects

{% assign featured_projects = site.projects | where: "featured", true | sort: "date" | reverse %}

{% if featured_projects.size > 0 %}
<div class="projects-grid">
  {% for project in featured_projects %}
  <div class="project-card">
    {% if project.thumbnail %}
    <div class="project-thumbnail">
      <img src="{{ project.thumbnail | relative_url }}" alt="{{ project.title }}">
    </div>
    {% endif %}
    <div class="project-info">
      <h3><a href="{{ project.url | relative_url }}">{{ project.title }}</a></h3>
      <p class="project-description">{{ project.description }}</p>
      {% if project.technologies %}
      <div class="project-tech">
        {% for tech in project.technologies limit:5 %}
          <span class="tech-badge">{{ tech }}</span>
        {% endfor %}
      </div>
      {% endif %}
      <p class="project-meta">
        {% if project.status %}<span class="status-{{ project.status | downcase | replace: ' ', '-' }}">{{ project.status }}</span>{% endif %}
        {% if project.date %} | {{ project.date | date: "%B %Y" }}{% endif %}
      </p>
      <a href="{{ project.url | relative_url }}" class="read-more">View Project →</a>
    </div>
  </div>
  {% endfor %}
</div>
{% else %}
<p><em>Featured projects coming soon...</em></p>
{% endif %}

## All Projects

{% assign all_projects = site.projects | sort: "date" | reverse %}

{% if all_projects.size > 0 %}
<div class="projects-list">
  {% for project in all_projects %}
  <div class="project-list-item">
    <h3><a href="{{ project.url | relative_url }}">{{ project.title }}</a></h3>
    <p class="project-date">{{ project.date | date: "%B %Y" }}</p>
    <p>{{ project.description }}</p>
    {% if project.domain %}
    <p class="project-domains">
      <strong>Domain:</strong> {{ project.domain | join: ', ' }}
    </p>
    {% endif %}
  </div>
  {% endfor %}
</div>
{% else %}
<p><em>Projects will be added soon. Check back later!</em></p>
{% endif %}

---

## Filter by Domain

*Coming soon: Filter projects by domain, technology, or status*

## Project Categories

- **Bioinformatics**: scRNA-seq, genomics, systems biology
- **Machine Learning**: Classification, clustering, deep learning
- **Physics & Scientific Computing**: Quantum systems, topological analysis
- **Data Engineering**: Pipelines, ETL, infrastructure
- **Visualization**: Interactive dashboards, scientific figures
