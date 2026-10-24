# T-OCS V7 — Public L0 Partner Overview

This folder contains the public-facing, non-confidential overview for T-OCS V7.0 REVISION-2.

- Site source: [index.html](./index.html)
- **Temporary public preview (not production hosting):** https://htmlpreview.github.io/?https://github.com/NAIF-Automotive-Safety-Lab/Car-seat-engineering/blob/main/docs/inbound/index.html
- Project contact: naifautomotivesafetylab@gmail.com
- Partner inquiry guide: [PARTNER_INQUIRY.md](./PARTNER_INQUIRY.md)
- Optional public intake: [GitHub non-confidential partner inquiry](https://github.com/NAIF-Automotive-Safety-Lab/Car-seat-engineering/issues/new?template=partner-capability.yml). GitHub issues are public; do not submit confidential or proprietary information.

## Disclosure and evidence boundary

NAIF Automotive Safety Lab is an independent automotive occupant-safety research and engineering initiative under SILVER TRUCKERS LLC. This page describes a rear-seat occupant-control system concept in engineering definition / evidence closure.

The physical test article, fixture, dynamic/sled testing, measurement/correlation and validation are future or partner-dependent steps unless separately established by evidence. This page does not claim validated safety performance, certification, regulatory compliance, production readiness, OEM endorsement, patent status, or customer traction. No proprietary dimensions, mechanism details, or confidential geometry should be added to this public page.

The contact form is a mailto-only interaction. It does not submit or store visitor data on the website; the visitor must review and send the draft in their own email application. No inquiry form has been submitted during implementation and no email has been sent.

## Latest verified QA — 24 October 2026

- The temporary preview returned HTTP 200 and rendered the intended T-OCS V7 page.
- English / Arabic toggle passed a round-trip test. Arabic produced `lang=ar`, `dir=rtl`, and the Arabic H1; English restoration produced `lang=en`, `dir=ltr`, and the English H1.
- No horizontal overflow was observed at 320px (body scrollWidth 300px, clientWidth 316px), 360px or 390px.
- The live Replit deployment is **not ready**: the status endpoint may report `success`, but an external fetch still returned a Replit build placeholder rather than the real site.
- GitHub Pages production URL still returns HTTP 404. The connected GitHub Actions integration lacked permission to enable/create the Pages site.

## Hosting status

The preview link above uses a third-party HTML preview service to display the public repository file. It is a temporary preview route, not a branded production-domain deployment or uptime guarantee. Do not share it as a permanent official domain. A production deployment should use the approved lab-controlled hosting/domain after the deployed page and its interactions are verified.