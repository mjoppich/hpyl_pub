/* global window document */

import * as React from 'react';
import { render } from 'react-dom';
import { BrowserRouter as Router } from 'react-router-dom';
import { MainApp } from './pages/App';
import config from './main_pages/config';


const AppClient = () => 
<Router basename={'/'+config.restFolder}>
  <MainApp/>
</Router>
;

window.onload = () => {
  render(<AppClient />, document.getElementById('main'));
};
