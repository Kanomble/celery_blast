from django.test import SimpleTestCase

from blast_project.result_responses import generated_html_response


class GeneratedResultResponseTests(SimpleTestCase):
    def test_static_result_response_uses_restrictive_headers(self):
        response = generated_html_response('<p>plain result</p>')

        self.assertEqual('text/html; charset=utf-8', response['Content-Type'])
        self.assertEqual('nosniff', response['X-Content-Type-Options'])
        self.assertEqual('DENY', response['X-Frame-Options'])
        self.assertIn("default-src 'none'", response['Content-Security-Policy'])
        self.assertIn('sandbox allow-downloads', response['Content-Security-Policy'])
        self.assertNotIn('allow-scripts', response['Content-Security-Policy'])

    def test_interactive_result_response_sandboxes_visualization_scripts(self):
        html = '<script>window.rendered = true;</script><div>visualization</div>'
        response = generated_html_response(html, interactive=True)

        self.assertContains(response, '<script>window.rendered = true;</script>', html=False)
        self.assertIn('sandbox allow-scripts', response['Content-Security-Policy'])
        self.assertNotIn('allow-same-origin', response['Content-Security-Policy'])
        self.assertIn('https://cdn.bokeh.org', response['Content-Security-Policy'])
