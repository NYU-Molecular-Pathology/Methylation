## Generate a REDCap API token
1. Login to REDCap and from the Projects homepage go to the "All Samples Database"
2. Next, click on the sidebar where it says "API" highlighted in red below.<br>
**Note**: *If you are missing this API access menu, access can be added for API  privilege under the "User Rights" menu in the same sidebar by someone who has access.  For new projects, access to the API module must be requested from the redcap admins.*<br>
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/00be73964a9542938a89810f074bb3ec03e724ed/screenshots/redcap_sidebar.png" alt="redcapsidebar" width="250"/><br />
3. On the following page below, click "Generate Token" and you can save this token locally only as it never expires and must be protected.
<img src="https://github.com/NYU-Molecular-Pathology/Methylation/blob/00be73964a9542938a89810f074bb3ec03e724ed/screenshots/redcap_api.png" alt="apidraw" width="850"/><br />
Note that if the REDCap version is updated the link may change from /redcap_v12.0.28/API to a different version.<br>
In this case, the project ID is 12352 so the url ends in pid=12352.  Different REDCap projects will a distinct REDCap token.
