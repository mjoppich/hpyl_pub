import * as React from 'react'; 

import ChipInput from 'material-ui-chip-input'
import AutoComplete from 'material-ui/AutoComplete'
import axios from 'axios';
import config from '../config';

import SelectField from 'material-ui/SelectField';
import MenuItem from 'material-ui/MenuItem';

export interface SwissModelViewerProps { alignments: any};
export interface SwissModelViewerState { selectedElement: any, resultPart:any};

export default class SwissModelViewer extends React.Component<SwissModelViewerProps, SwissModelViewerState>{

    neo4jd3: any = null;

    constructor(props)
    {
        super(props);
    }

    componentWillMount()
    {
        var self=this;
        this.setState({selectedElement: undefined, resultPart: ""});

    }

    componentDidMount(){

    }

    handleChange(event, index, value)
    {
        var newValue = value;

        var self = this;

        axios.post(config.getRestAddress() + "/swissmodel/query", {'uniprot': newValue}, config.axiosConfig)
        .then(function (response) {
            console.log(response.data)

            var newResultPart = self.makeResult(response.data);

            self.setState({resultPart:newResultPart})

        })
        .catch(function (error) {
          console.log(error)

          var newResultPart = <p>There is no SwissModel Information for protein: {newValue}</p>;

          self.setState({resultPart: newResultPart})
        });

        this.setState({selectedElement: newValue})
    }

    makeResult(inputData)
    {
        console.log("make result");
        console.log(inputData);

        var errorMessages = [];
        var structStuff = [];

        var uniprotIDs = Object.keys(inputData);

        for (var u = 0; u < uniprotIDs.length; ++u)
        {
            var queryID = uniprotIDs[u];
            var resultElem = inputData[queryID];

            if ((resultElem.result == undefined) || (resultElem.result.structures == undefined) || (resultElem.result.structures.length == 0))
            {
                errorMessages.push(<p key={errorMessages.length}>We tried finding one, but there is no SwissModel Structures for protein: {queryID}</p>);
            } else {
    
                var allStructures = resultElem.result.structures;
    
    
                for (var i = 0; i < allStructures.length; ++i)
                {
                    var strucElem = allStructures[i];
    
                    var modAlignment = strucElem.alignment.replace(/\\n/g, "\n");
                    var templateArr = strucElem.template.split('\.');
                    var templateID = strucElem.template;
                    
                    if (templateArr.length >= 2)
                    {
                        templateID = templateArr[0] + "." + templateArr[1];
                    }
    
                    var identity = parseFloat((strucElem.identity).toFixed(5));
                    var similarity = parseFloat((strucElem.similarity).toFixed(5));
                    var coverage = parseFloat((strucElem.coverage).toFixed(5));
                    var qMean = parseFloat((strucElem.qmean).toFixed(5));
                    var qmean_norm = parseFloat((strucElem.qmean_norm).toFixed(5));
                    var gmqe = parseFloat((strucElem.gmqe).toFixed(5));
    
                    var rowStyle = {};
                    if (structStuff.length % 2 === 0)
                    {
                        rowStyle = {
                            backgroundColor: config.table_even_bg
                        }
                    }
    
                    structStuff.push(
                        <tr key={structStuff.length} style={rowStyle}>
                            <td style={{verticalAlign: 'top'}}>{queryID}<br/>{strucElem.method}<br/>{strucElem['oligo-state']}<br/>{strucElem.provider}</td>
                            <td style={{verticalAlign: 'top'}}>{strucElem.from}-{strucElem.to}<br/>{identity}<br/>{similarity}<br/>{coverage}</td>
                            <td style={{verticalAlign: 'top'}}>{qMean}<br/>{qmean_norm}<br/>{gmqe}</td>
                            <td style={{verticalAlign: 'top', display: 'block', overflow:'scroll', width: '600px'}}><pre>{modAlignment}</pre></td>
                            <td style={{verticalAlign: 'top'}}><a target="_blank" href={"https://swissmodel.expasy.org/templates/"+templateID}>{strucElem.template}</a></td>
                        </tr>
                    );
                }
            }
        }

        var tableElems = <div></div>;

        if (structStuff.length > 0)
        {
            tableElems = <table style={{ tableLayout: "unset", width: "100%"}}>
                            <tbody>
                                <tr style={{textAlign: 'left'}}>
                                    <th>Method<br/>Oligo-State<br/>provider</th>
                                    <th>From/To<br/>Identity<br/>Similarity<br/>Coverage</th>
                                    <th>QMean/QMean-Norm/GMQE</th>
                                    <th>Alignment</th>
                                    <th>Template<br/>SMTL ID</th>
                                </tr>
                                {structStuff}
                            </tbody>
                        </table>;
        }

        return (
            <div>
            {tableElems}
            {errorMessages}
            </div>
        );


    }
  
    render(){

        console.log("SwissModelViewer")
        console.log(this.props.alignments);

        var menuitems = [];

        for (var i = 0; i < this.props.alignments.msa.length; ++i)
        {
            var rowAlign = this.props.alignments.msa[i];

            if ((rowAlign['xrefs']['Uniprot'] == null) || (rowAlign['xrefs']['Uniprot'] == undefined) ||(rowAlign['xrefs']['Uniprot'].length == 0))
            {
                continue;
            }

            var labelText = rowAlign['xrefs']['Uniprot'] + " (" + rowAlign.entryID + ")";
            console.log(labelText);

            var elemitem = <MenuItem key={i} value={rowAlign['xrefs']['Uniprot'][0]} primaryText={labelText} />
            menuitems.push(elemitem);
        }

        var menuPart =  <SelectField hintText="Select one or multiple proteins for SwissModel lookup." multiple={true} value={this.state.selectedElement} onChange={this.handleChange.bind(this)}>
                            {menuitems}
                        </SelectField>

        return (<div>
            
            <h3>Select the protein:</h3>
            {menuPart}
            {this.state.resultPart}

            </div>);
    }
}