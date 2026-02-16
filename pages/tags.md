---
layout: page
title: "Posts by Tag"
description: "Browse blog posts organized by tag"
permalink: /blog/tags/
---

# Posts by Tag

<div class="tag-cloud">
{% assign tags = site.tags | sort %}
{% for tag in tags %}
  <a href="#{{ tag[0] | slugify }}" class="tag-cloud-item">{{ tag[0] }} ({{ tag[1].size }})</a>
{% endfor %}
</div>

---

{% for tag in tags %}
<div class="tag-section" id="{{ tag[0] | slugify }}">

## {{ tag[0] }}

<ul class="tag-posts">
{% for post in tag[1] %}
  <li>
    <a href="{{ post.url | relative_url }}">{{ post.title }}</a>
    <span class="post-date">{{ post.date | date: "%B %d, %Y" }}</span>
  </li>
{% endfor %}
</ul>

</div>
{% endfor %}

{% if site.tags.size == 0 %}
<p><em>No tags yet. Tags will appear here as blog posts are published.</em></p>
{% endif %}
